"""Simple tool for displaying DB connections, used for debugging performance issues
"""
from __future__ import print_function, division
import os, sys, subprocess, re
from collections import OrderedDict
from datetime import tzinfo, timedelta, datetime
from multipatch_analysis import config
from sqlalchemy import create_engine


# https://stackoverflow.com/questions/796008/cant-subtract-offset-naive-and-offset-aware-datetimes
class UTC(tzinfo):
  def utcoffset(self, dt):
    return timedelta(0)
  def tzname(self, dt):
    return "UTC"
  def dst(self, dt):
    return timedelta(0)

utc = UTC()


engine = create_engine(config.synphys_db_host_rw + '/postgres?application_name=db_monitor')
conn = engine.connect()


def list_db_connections():
    now = datetime.now(utc)
    known_addrs = config.known_addrs
    
    tr = conn.begin()
    query = conn.execute("select client_addr, client_port, datname, query, state, usename, application_name, pid, state_change from pg_stat_activity;")
    result = query.fetchall()
    tr.rollback()
    
    recs = []
    for r in result:
        r = OrderedDict(r.items())
        ip = r['client_addr']
        r['hostname'] = None if ip is None else hostname(ip)
        r['user'] = known_addrs.get(ip, known_addrs.get(r['hostname'], None))
        recs.append(r)
        
        last_change = r['state_change']
        if last_change is None:
            r['idle_time'] = None
        else:
            r['idle_time'] = (now - last_change).total_seconds()
    
    return recs


def kill_pid(pid, terminate=False):
    tr = conn.begin()
    if terminate:
        result = conn.execute("select pg_terminate_backend(%d);" % pid).fetchall()
    else:
        result = conn.execute("select pg_cancel_backend(%d);" % pid).fetchall()
    tr.rollback()
    return result[0][0]


_known_hostnames = {}
def hostname(ip):
    global _known_hostnames
    if ip in _known_hostnames:
        return _known_hostnames[ip]
    try:
        host = subprocess.check_output(['host', ip]).partition('pointer ')[2].rstrip('.\n')
    except subprocess.CalledProcessError:
        host = "hostname not found"
    _known_hostnames[ip] = host
    return host


def check():
    print("====================  DB connections ======================")    
    connects = list_db_connections()
    connect_users = [(conn['client_addr'], conn['hostname'], conn['user']) for conn in connects]
    
    users = list(set(connect_users))
    counts = {user:connect_users.count(user) for user in users}
    users.sort(key=lambda user: counts[user], reverse=True)
    
    for user in users:
        (ip, host, name) = user
        host = host or "[None]"
        name = name or '???'
        count = counts[user]
        print("{:10s}{:15s}".format(name, ip))
        
        for con in connects:
            if con['client_addr'] != ip:
                continue
            app = con['application_name']
            for pkg in ['acq4', 'multipatch_analysis']:
                a,b,c = app.partition(pkg)
                if b == '':
                    continue
                else:
                    app = c
                    break
                    
            if con['idle_time'] is None:
                state = '[%s]' % con['state']
            else:
                state = '[%s %0.1fhr]' % (con['state'], con['idle_time'] / 3600.)
            
            query = con['query'].replace('\n', ' ')[:110]
            
            print("          {:15s} {:15s} {:45s} {:6d} {:16s} {:s}   ".format(con['usename'], con['datname'], app[:45], con['pid'], state, query))
    

class ScreenPrint(object):
    """Simple interface for printing full-screen text
    """
    def __init__(self):
        global print
        self._print = print
        print = self.print_line
        
        self.current_row = 0
        self.rows = 100
        self.columns = 100

    def reset(self):
        self.current_row = 0
        self._print("\033[0;0f")
        self.rows, self.columns = list(map(int, os.popen('stty size', 'r').read().split()))

    def print_line(self, msg):
        while len(msg) > self.columns:
            line, msg = msg[:self.columns], msg[self.columns:]
            self._print(line)
            self.current_row += 1
        self._print(("{:%ds}" % self.columns).format(msg))
        self.current_row += 1
        
    def clear_to_bottom(self):
        for i in range(self.current_row, self.rows-1):
            self.print_line("")


if __name__ == '__main__':
    # import pyqtgraph as pg
    # pg.dbg()
    import time, argparse
    
    parser = argparse.ArgumentParser(description="Monitor and manipulate postgres database connections")
    parser.add_argument('--one', action='store_true', default=False, help="List all connections once, then exit.")
    parser.add_argument('--kill', action='store_true', default=False, help="Close idle connections")
    parser.add_argument('--terminate', action='store_true', default=False, help="Terminate (rather than cancel) connections")
    parser.add_argument('--kill-from', type=str, default=None, help="Kill only connections from a particular user or IP address", dest='kill_from')
    parser.add_argument('--kill-age', type=float, default=0.5, help="Kill only connections idle for some number of hours (default=0.5)", dest='kill_age')
    args = parser.parse_args()
    
    if args.kill:
        to_kill = list_db_connections()
        to_kill = [c for c in to_kill if c['idle_time']/3600. > args.kill_age]
        if args.kill_from is not None:
            to_kill = [c for c in to_kill if args.kill_from in (c['user'], c['client_addr'], c['hostname'])]
        for c in to_kill:
            killed = kill_pid(c['pid'], args.terminate)
            print("%d %s" % (c['pid'], 'killed' if killed else 'kill failed'))
    
    if args.one:
        check()
    else:
        sp = ScreenPrint()
        while True:
            sp.reset()
            check()
            sp.clear_to_bottom()        
            time.sleep(1.0)


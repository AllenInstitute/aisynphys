import os, sys, subprocess, time, re

argv = sys.argv[1:]

usage = """
Run pair_analysis.py on a list of pair IDs
Usage:  python pair_analysis_list.py expt_id1 pre_cell_id1 post_cell_id1 ...
or      python pair_analysis_list.py - 
 (in the second case, send pair IDs via stdin, one ID per line)
"""

if len(argv) == 0:
    print(usage)
    sys.exit(0)

def get_next_pair():
    global argv
    if len(argv) == 0:
        return None
    if argv[0] == '-':
        while True:
            line = sys.stdin.readline()
            if line == '':
                return None
            if '#' in line:
                line = line.partition('#')[0]
            line = line.strip()
            if len(line) == 0:
                continue
            return re.split(r'\s+', line)
    else:
        return (argv.pop(0), argv.pop(0), argv.pop(0))


procs = []
pa_script = os.path.join(os.path.dirname(__file__), 'pair_analysis.py')


while True:
    if len(procs) < 2:
        syn_id = get_next_pair()
        if syn_id is not None:
            proc = subprocess.Popen(['python', pa_script, syn_id[0], syn_id[1], syn_id[2]])
            procs.append(proc)
    # remove finished processes from list
    procs = [p for p in procs if p.poll() is None]
    time.sleep(0.1)
    
    if len(procs) == 0:
        break

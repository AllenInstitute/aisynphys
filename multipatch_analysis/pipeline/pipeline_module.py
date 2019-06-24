from __future__ import division, print_function
import sys, time, multiprocessing, traceback
from datetime import datetime
import numpy as np
from collections import OrderedDict
from pyqtgraph import toposort
from .. import database as db


class PipelineModule(object):
    """Pipeline modules represent analysis tasks that can be run independently of other parts of the analysis
    pipeline. 
    
    For any given experiment, a sequence of analysis stages must be processed in order. Each stage requires a specific
    set of inputs to be present, and produces/stores some output. Inputs and outputs can be raw data files, database tables,
    etc., although outputs are probably always written into a database. Each PipelineModule subclass represents a single stage
    in the sequence of analyses that occur across a pipeline.
    
    The work done by a single stage in the pipeline is divided up into jobs (units of work), where each stage may
    decide for itself what a suitable unit of work is. For many stages, the unit of work is the analysis done on a single experiment.
    Some stages may have finer granularity (multiple jobs per experiment), and some stages may have coarser granularity
    (multiple experiments per job), or even only have a single unit of work for the entire database, such as when aggregating 
    very high-level results.
    
    Note that PipelineModule classes generally need not implement any actual _analysis_; rather, they are simply responsible for
    data _management_ within the analysis pipeline. When an PipelineModule is asked to update an analysis result, it may call out to
    other packages to do the real work.
    
    From here, we should be able to:
    - Ask for which jobs this analysis has already been run (and when)
    - Ask for which jobs this analysis needs to be run (input dependencies are met and no result exists yet or is out of date)    
    - Run and store analysis for any specific experiment
    - Remove / re-run analysis for any specific experiment
    - Remove / re-run analysis for all jobs (for example, after a code change affecting all analysis results)
    - Initialize storage for analysis results (create DB tables, empty folders, etc.)
    
    """    
    
    name = None
    dependencies = []
    maxtasksperchild = None

    @staticmethod
    def all_modules():
        """Return an ordered dictionary mapping {name: module} of all known pipeline modules,
        topologically sorted by dependencies (least dependent to most dependent).
        """
        subclasses = [PipelineModule]
        for sc in subclasses:
            subclasses.extend(sc.__subclasses__())
        
        excluded = [PipelineModule, DatabasePipelineModule]
        deps = {c:c.dependencies for c in subclasses if c not in excluded}
        return OrderedDict([(mod.name, mod) for mod in toposort(deps)])
    
    @classmethod
    def dependent_modules(cls):
        """Return a list of other modules that directly depend on this module.
        """
        mods = cls.all_modules()
        return [mod for mod in mods.values() if cls in mod.dependencies]
    
    @classmethod
    def all_dependent_modules(cls):
        """Return a list of other modules that directly or indirectly depend on this module.
        """
        deps = cls.dependent_modules()
        for dep in deps:
            deps.extend(dep.dependent_modules())
        return [mod for mod in PipelineModule.all_modules().values() if mod in deps]
    
    @classmethod
    def update(cls, job_ids=None, retry_errors=False, limit=None, parallel=False, workers=None, raise_exceptions=False):
        """Update analysis results for this module.
        
        Parameters
        ----------
        job_ids : list | None
            List of job IDs to be updated, or None to update all jobs.
        retry_errors : bool
            If True, jobs that previously failed will be attempted again.
        parallel : bool
            If True, run jobs in parallel threads or subprocesses.
        workers : int or None
            Number of parallel workers to spawn. If None, then use one worker per CPU core.
        limit : int | None
            Maximum number of jobs to process (or None to disable this limit).
            If limit is enabled, then jobs are randomly shuffled before selecting the limited subset.
        raise_exceptions : bool
            If True, then exceptions are raised and will end any further processing.
            If False, then errors are logged and ignored.
            This is used mainly for debugging to allow traceback inspection.
        """
        print("Updating pipeline stage: %s" % cls.name)
        n_retry = 0
        if job_ids is None:
            print("Searching for jobs to update..")
            drop_job_ids, run_job_ids, error_job_ids = cls.updatable_jobs()
            
            if retry_errors:
                run_job_ids += error_job_ids
                n_retry = len(error_job_ids)
            
            if limit is not None:
                # pick a random subset to import; this is just meant to ensure we get a variety
                # of data when testing the import system.
                rng = np.random.RandomState(0)
                rng.shuffle(run_job_ids)
                run_job_ids = run_job_ids[:limit]
                drop_job_ids = [jid for jid in drop_job_ids if jid in run_job_ids]
        else:
            run_job_ids = job_ids
            drop_job_ids = job_ids
            
        print("Found %d jobs to update." % len(run_job_ids))

        # drop invalid records first
        if len(drop_job_ids) > 0:
            print("Dropping %d invalid results (will not update).." % len(drop_job_ids))
            print(drop_job_ids)
            cls.drop_jobs(drop_job_ids)
        if len(run_job_ids) > 0:
            print("Dropping %d invalid results (will update).." % len(run_job_ids))
            cls.drop_jobs(run_job_ids)

        run_jobs = [(expt_id, i, len(run_job_ids)) for i, expt_id in enumerate(run_job_ids)]

        if parallel:
            # kill DB connections before forking multiple processes
            db.dispose_engines()
            
            print("Processing all jobs (parallel)..")
            pool = multiprocessing.Pool(processes=workers, maxtasksperchild=cls.maxtasksperchild)
            try:
                # would like to just call cls._run_job, but we can't pass a method to Pool.map()
                # instead we wrap this with the run_job_parallel function defined below.
                parallel_jobs = [(cls, job) for job in run_jobs]
                job_results = {}
                chunksize = cls.maxtasksperchild or 1
                for result in pool.imap(run_job_parallel, parallel_jobs, chunksize=chunksize):  # note: maxtasksperchild is broken unless we also force chunksize
                    job_results[result['job_id']] = result['error']
                    print("Finished %d/%d  (%0.1f%%)" % (len(job_results), len(run_jobs), 100*len(job_results)/len(run_jobs)))
            finally:
                pool.close()
                
        else:
            print("Processing all jobs (serial)..")
            job_results = {}
            for job in run_jobs:
                result = cls._run_job(job, raise_exceptions=raise_exceptions)
                job_results[result['job_id']] = result['error']
                
        errors = {job:result for job,result in job_results.items() if result is not None}
        return {'n_dropped': len(drop_job_ids), 'n_updated': len(run_job_ids), 'n_errors': len(errors), 'errors': errors, 'n_retry': n_retry}

    @classmethod
    def _run_job(cls, job, raise_exceptions=False):
        """Entry point for running a single analysis job; may be invoked in a subprocess.
        """
        job_id, job_index, n_jobs = job
        print("Processing %s %d/%d  %0.3f" % (cls.name, job_index+1, n_jobs, job_id))
        start = time.time()
        try:
            cls.process_job(job_id)
        except Exception as exc:
            if raise_exceptions:
                raise
            else:
                print("Error processing %s %d/%d  %0.3f:" % (cls.name, job_index+1, n_jobs, job_id))
                sys.excepthook(*sys.exc_info())
                return {'job_id': job_id, 'error': str(exc)}
        else:
            print("Finished %s %d/%d  %0.3f  (%0.2f sec)" % (cls.name, job_index+1, n_jobs, job_id, time.time()-start))
            return {'job_id': job_id, 'error': None}
    
    @classmethod
    def process_job(cls, job_id):
        """Process analysis for one job.
        
        Parameters
        ----------
        job_id :
            The ID of the job to be processed.
        
        Must be reimplemented in subclasses.
        """
        raise NotImplementedError()

    @classmethod
    def initialize(cls):
        """Create space (folders, tables, etc.) for this analyzer to store its results.
        """
        raise NotImplementedError()
        
    @classmethod
    def drop_jobs(cls, job_ids):
        """Remove all results previously stored for a list of job IDs.
        """
        raise NotImplementedError()

    @classmethod
    def drop_all(cls):
        """Remove all results generated by this module.
        """
        raise NotImplementedError()

    @classmethod
    def finished_jobs(cls):
        """Return an ordered dict of job IDs that have been successfully processed by this module,
        the dates when they were processed, and whether each job succeeded:  {job_id: (date, success)}

        Note that some results returned may be obsolete if dependencies have changed.
        """
        raise NotImplementedError()

    @classmethod
    def ready_jobs(cls):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        # default implpementation collects IDs of finished jobs from upstream modules.
        job_times = OrderedDict()
        for i,mod in enumerate(cls.dependencies):
            jobs = mod.finished_jobs()
            for job_id,(ts,success) in jobs.items():
                if success is False:
                    continue
                job_times.setdefault(job_id, [None]*len(cls.dependencies))
                job_times[job_id][i] = ts
            
        ready = OrderedDict()
        for job_id, times in job_times.items():
            if None in times:
                continue
            ready[job_id] = max(times)
            
        return ready

    @classmethod
    def updatable_jobs(cls):
        """Return lists of jobs that should be updated and/or should have their results dropped.
        
        Returns
        -------
        drop_job_ids : list
            Jobs that need to be dropped because their output is invalid and they are not ready to be updated again
        run_job_ids : list
            Jobs that need to be updated AND are ready to be updated
        error_job_ids : list
            Jobs that previously failed to update due to an error
        """
        run_job_ids = []
        drop_job_ids = []
        error_job_ids = []
        ready = cls.ready_jobs()
        finished = cls.finished_jobs()
        for job in ready:
            if job in finished:
                date, success = finished[job]
                if ready[job] > date:
                    # result is invalid
                    run_job_ids.append(job)
                else:
                    if success is False:
                        error_job_ids.append(job)
                    else:
                        # result is valid
                        pass  
            else:
                # no current result
                run_job_ids.append(job)
        
        # look for orphaned results
        for job in finished:
            if job not in ready and job not in drop_job_ids:
                drop_job_ids.append(job)
        
        print("%d jobs ready for processing, %d finished, %d need drop, %d need update, %d previous errors" % (len(ready), len(finished), len(drop_job_ids), len(run_job_ids), len(error_job_ids)))
        return drop_job_ids, run_job_ids, error_job_ids


def run_job_parallel(job):
    # multiprocessing Pool.map doesn't work on methods; must be a plain function
    cls, job = job
    return cls._run_job(job)


class DatabasePipelineModule(PipelineModule):
    """PipelineModule that implements default behaviors for interacting with database.
    
    * Manages adding / removing entries from pipeline table to indicate job status
    * Default implementations for dropping records 
    """
    table_group = None

    @classmethod
    def create_db_entries(cls, job_id, session):
        """Generate DB entries for *job_id* and add them to *session*.
        """
        raise NotImplementedError()
        
    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        raise NotImplementedError()

    @classmethod
    def dependent_job_ids(cls, module, job_ids):
        """Return a list of all finished job IDs in this module that depend on 
        specific jobs from another module.
        """
        if module not in cls.dependencies:
            raise ValueError("%s does not depend on module %s" % (cls, module))
        
        # In most cases, modules use the same IDs as the modules that they depend on.
        # (usually this is the experiment ID)
        return job_ids

    @classmethod
    def process_job(cls, job_id):
        session = db.Session(readonly=False)
        # drop old pipeline job record
        session.query(db.Pipeline).filter(db.Pipeline.job_id==job_id).filter(db.Pipeline.module_name==cls.name).delete()
        session.commit()
        
        try:
            errors = cls.create_db_entries(job_id, session)
            job_result = db.Pipeline(module_name=cls.name, job_id=job_id, success=True, error=errors, finish_time=datetime.now())
            session.add(job_result)

            session.commit()
        except Exception:
            session.rollback()
            
            err = ''.join(traceback.format_exception(*sys.exc_info()))
            job_result = db.Pipeline(module_name=cls.name, job_id=job_id, success=False, error=err, finish_time=datetime.now())
            session.add(job_result)
            session.commit()
            raise
        finally:
            session.close()

    @classmethod
    def initialize(cls):
        """Create space (folders, tables, etc.) for this analyzer to store its results.
        """
        print("Initialize tables in %s" % cls.name)
        cls.table_group.create_tables()
        
    @classmethod
    def drop_all(cls, reinitialize=True, drop_deps=True):
        """Remove all results generated by this module and dependent modules.
        """
        if drop_deps:
            all_mods = cls.all_dependent_modules()
            for dep in reversed(all_mods):
                dep.drop_all(reinitialize=False, drop_deps=False)
            
        # drop tables and pipeline job records
        print("Drop all in %s" % cls.name)
        cls.table_group.drop_tables()
        session = db.Session(readonly=False)
        session.query(db.Pipeline).filter(db.Pipeline.module_name==cls.name).delete()
        session.commit()

        if reinitialize:
            cls.initialize()
            for dep in all_mods:
                dep.initialize()
    
    @classmethod
    def drop_jobs(cls, job_ids, session=None, skip=None):
        """Remove all results previously stored for a list of job IDs.
        
        The associated results of dependent modules are also removed.
        """
        if session is None:
            session = db.Session(readonly=False)
        
        if skip is None:
            skip = []
        
        for dep in reversed(cls.dependent_modules()):
            if dep in skip:
                continue            
            dep_jobs = dep.dependent_job_ids(cls, job_ids)
            dep.drop_jobs(dep_jobs, session=session, skip=skip)
        
        print("Dropping %d jobs from %s module.." % (len(job_ids), cls.name))
        records = cls.job_records(job_ids, session)
        if len(records) == 0:
            print("   (no records to remove for these job IDs)")
        else:
            for i,rec in enumerate(records):
                session.delete(rec)
                print("   record %d/%d\r" % (i, len(records)), end='')
                sys.stdout.flush()
            session.query(db.Pipeline).filter(db.Pipeline.module_name==cls.name).filter(db.Pipeline.job_id.in_(job_ids)).delete(synchronize_session=False)
            print("   dropped %d records; committing.." % len(records))
            session.commit()
        
        skip.append(cls)  # only process each module once
    
    @classmethod
    def finished_jobs(cls):
        """Return an ordered dict of job IDs that have been successfully processed by this module,
        the dates when they were processed, and whether each job succeeded:  {job_id: (date, success)}

        Note that some results returned may be obsolete if dependencies have changed.
        """
        session = db.Session()
        jobs = session.query(db.Pipeline.job_id, db.Pipeline.finish_time, db.Pipeline.success).filter(db.Pipeline.module_name==cls.name).all()
        session.rollback()
        return OrderedDict([(uid, (date, success)) for uid, date, success in jobs])

"""
Script used for collecting stochastic release model results and running dimensionality reduction on them.

Depending on the number of synapses to process and the size of the model parameter space, this
may consume a large amount of memory and CPU time.

"""

import os, sys, pickle, gc, time, traceback, threading
import numpy as np
import umap
import sklearn.preprocessing, sklearn.decomposition
from aisynphys.stochastic_release_model import StochasticModelRunner, load_cached_model_results
from aisynphys.database import default_db as db
from aisynphys import config
from aisynphys.ui.progressbar import ProgressBar


class ThreadTrace(object):
    """ 
    Used to debug freezing by starting a new thread that reports on the 
    location of other threads periodically.
    """
    def __init__(self, interval=10.0):
        self.interval = interval
        self._stop = False
        self.start()

    def start(self, interval=None):
        if interval is not None:
            self.interval = interval
        self._stop = False
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True
        self.thread.start()

    def run(self):
        while True:
            if self._stop is True:
                return
                    
            print("\n=============  THREAD FRAMES:  ================")
            for id, frame in sys._current_frames().items():
                if id == threading.current_thread().ident:
                    continue

                # try to determine a thread name
                try:
                    name = threading._active.get(id, None)
                except:
                    name = None
                if name is None:
                    try:
                        # QThread._names must be manually set by thread creators.
                        name = QtCore.QThread._names.get(id)
                    except:
                        name = None
                if name is None:
                    name = "???"

                print("<< thread %d \"%s\" >>" % (id, name))
                traceback.print_stack(frame)
            print("===============================================\n")
            
            time.sleep(self.interval)

# tt = ThreadTrace(60)
print("start thread trace..")

cache_path = config.cache_path
model_result_cache_path = os.path.join(cache_path, 'stochastic_model_results')

cache_files = [os.path.join(model_result_cache_path, f) for f in os.listdir(model_result_cache_path)]

for run_type in ['all_results', 'likelihood_only']:
    ## Load all model outputs into a single array
    agg_result, cache_files, param_space = load_cached_model_results(cache_files, db=db)
    agg_shape = agg_result.shape
    print("  cache loaded.")

    if run_type == 'all_results':
        # first time, use complete result array
        flat_result = agg_result.reshape(agg_shape[0], np.product(agg_shape[1:]))
    elif run_type == 'likelihood_only':
        # next use only mode likelihood; no amplitude
        flat_result = agg_result[...,0].reshape(agg_shape[0], np.product(agg_shape[1:-1]))
    else:
        raise ValueError(run_type)
    print("  reshape.")

    # Prescale model data
    print("   Fitting prescaler...")
    scaler = sklearn.preprocessing.StandardScaler()
    n_obs = flat_result.shape[0]
    chunk_size = 10
    n_chunks = n_obs // chunk_size

    scaler.fit(flat_result)
    print("   Prescaler transform...")
    scaled = scaler.transform(flat_result)

    print("   Prescaler done.")

    # free up some memory
    del agg_result
    del flat_result
    gc.collect()

    print("free memory")

# fit standard PCA   (uses ~2x memory of input data)
#try:
#    start = time.time()
#    print("Fitting PCA...")
#    n_pca_components = 30#500
#    pca = sklearn.decomposition.PCA(n_components=n_pca_components)
#    pca.fit(scaled)
#    print("  PCA fit complete.")
#
#    # run PCA
#    print("PCA transform...")
#    pca_result = pca.transform(scaled)
#    pca_file = os.path.join(cache_path, 'pca.pkl')
#    pickle.dump({'result': pca_result, 'params': param_space, 'cache_files': cache_files, 'pca': pca}, open(pca_file, 'wb'))
#    print("   PCA transform complete: %s" % pca_file)
#except Exception as exc:
#    print("PCA failed:")
#    traceback.print_exc()
#finally:
#    print("PCA time: %d sec" % int(time.time()-start))


# umap  (uses ~1x memory of input data)
#try:
#    start = time.time()
#    n_umap_components = 15#32
#    reducer = umap.UMAP(
#        n_components=n_umap_components,
#    #     n_neighbors=5,
#        low_memory=False,
#        init='spectral',   # also try 'random'
#        verbose=True,
#    )
#
#    print("Fit UMAP...")
#    reducer.fit(scaled)
#    print("   UMAP fit done.")
#
#    print("UMAP transform...")
#    umap_result = reducer.transform(scaled)
#    umap_file = os.path.join(cache_path, 'umap.pkl')
#    pickle.dump({'result': umap_result, 'params': param_space, 'cache_files': cache_files}, open(umap_file, 'wb'))
#    pickle.dump(reducer, open('umap_model.pkl', 'wb'))
#    print("   UMAP transform complete: %s" % umap_file)
#except Exception as exc:
#    print("UMAP failed:")
#    traceback.print_exc()
#finally:
#    print("UMAP time: %d sec" % int(time.time()-start))



    # fit sparse PCA  (uses ~6x memory of input data)
    try:
        start = time.time()
        print("Fitting sparse PCA...")
        n_pca_components = 50
        pca = sklearn.decomposition.MiniBatchSparsePCA(n_components=n_pca_components, n_jobs=-1)
        pca.fit(scaled)
        print("  Sparse PCA fit complete.")
   
        # run sparse PCA
        print("Sparse PCA transform...")
        sparse_pca_result = pca.transform(scaled)
        sparse_pca_file = os.path.join(cache_path, f'sparse_pca_{run_type}.pkl')
        pickle.dump({'result': sparse_pca_result, 'params': param_space, 'cache_files': cache_files, 'sparse_pca': pca}, open(sparse_pca_file, 'wb'))
        print("   Sparse PCA transform complete: %s" % sparse_pca_file)
    except Exception as exc:
        print("Sparse PCA failed:")
        traceback.print_exc()
    finally:
        print("Sparse PCA time: %d sec" % int(time.time()-start))


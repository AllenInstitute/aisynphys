# coding: utf8
from __future__ import print_function, division

"""
So here's an interesting hypothesis: some cortical neurons "prefer" to form
bidirectional synaptic connections. That is, cells will connect together randomly,
but as soon as a connection A->B forms, there is a higher probability that the
reciprocal connection B->A will form.

How could we test this hypothesis? Measure the actual fraction of connections
between many pairs of neurons, and assume this is a good estimate of the underlying
connection probability (Pc), then further assume that under the null hypothesis
(cells do have no preference for making or avoiding reciprocal connections), the
underlying reciprocal connection probability (Pr) is just Pr = Pc^2. Now we measure
the fraction of reciprocal connections found in our sample, and compare that 
to the estimated Pr. 

There are two major problems with this approach:

1. The simple relationship Pr = Pc^2 is only true in general if the probability
   of connection between all cell pairs sampled is uniform. If we are actually
   sampling from multiple cell types with unequal connection probabilities, then
   Pc^2 may be a bad estimate and our null hypothesis is invalid.
   (see https://www.mitpressjournals.org/doi/10.1162/NETN_a_00004)
   
2. Statistical power is very expensive here--if we want to state with good
   confidence that our measured reciprocal connectivity is different from the
   estimated Pr, then we require a very large N.


The tools below are meant to simulate a typical series of multipatch experiments
such that we can explore how the underlying cell type-dependent connection
probabilities (Wab) affect the measured reciprocal connectivity, and how 
sample size affects statistical confidence.



"""

import numpy as np


def run_expt(Wab, n_cells=4, n_expts=100, n_trials=1000):
    """Simulate many trials of a series of multipatch experiments.
    
    Parameters
    ----------
    Wab : array
        An NxN matrix specifying the probability of a connection forming between
        different types of cells. The size N of the matrix is the number of
        cell types, and each element [i,j] gives the probability (0.0-1.0) that
        a cell of type i will form a connection onto a cell of type j.
    n_cells : int
        The number of cells to test in each simulated multipatch experiment
        (default is 4). For each experiment, we select this many cells randomly
        from the N cell types. Connections are assigned randomly between cells
        according to the probabilities in Wab.
    n_experiments : int
        The number of simulated multipatch experiments to run before the series
        is complete. After this many experiments, the experimenter would gether
        their data and attempt to decide whether the reciprocal connectivity in
        their sample is higher or lower than random chance.
    n_trials : int
        The number of times to repeat the series of n_experiments. This allows
        us to evaluate the reproducibility of the results.
        
    Returns
    -------
    results : array
        An array of all experiment results. The shape is (n_trials, n_expts),
        and each element contains the results of a single experiment.
        the dtype is compound with three fields:
        
            * conn : the number of connections found
            * recip : the number of reciprocal connections found
            * probed : the number of connections probed
            
        For example, the total connection probability for the 5th trial can be
        calculated as ``results[4]['conn'].sum() / results[4]['probed'].sum()``.
    """
    Wab = np.array(Wab)
    print("-----------")
    print("cells: %d   expts: %d   trials: %d" % (n_cells, n_expts, n_trials))
    print(Wab)
    n_cell_types = Wab.shape[0]
    results = np.empty((n_trials, n_expts), dtype=[('conn', int), ('recip', int), ('probed', int)])
    
    # Run many trials so we can analyze the trial-to-trial variability in results
    for i in range(n_trials):
        # Run many experiments to recover the connection probability
        # (this is like running many multipatch experiments on a single class of neuron)
        for j in range(n_expts):
            # Randomly choose N cells from available cell types
            types = np.random.randint(n_cell_types, size=n_cells)
            
            # i,j connection probability matrix for this experiment
            cpm = Wab[types][:, types]
            
            # i,j boolean connection matrix for this experiment
            conn = cpm > np.random.random(size=(n_cells, n_cells))
            
            # clear diagonal
            diag = np.eye(n_cells, dtype='bool')
            conn[diag] = False
            
            # count total connections and reciprocal connections
            results[i, j]['conn'] = conn.sum()
            results[i, j]['recip'] = (conn & conn.T).sum()
            results[i, j]['probed'] = n_cells * (n_cells-1)

    cprobs = results['conn'].sum(axis=1) / results['probed'].sum(axis=1)
    rprobs = results['recip'].sum(axis=1) / results['probed'].sum(axis=1)
    ex_rprobs = cprobs**2
    ratios = rprobs / ex_rprobs
    
    # Connection probability we expect to measure, given Wab:
    print("  expected connection probability: %0.04f" % Wab.mean())
    # Actual connection probability measured across all trials:
    print("  measured connection probability: %0.4f±%0.4f" % (np.mean(cprobs), np.std(cprobs)))
    
    # reciprocal connection probability we expect to measure, given Wab:
    # print("  expected reciprocal probability: %0.4f" % (Wab**2).mean())  # this is not right
    # Actual reciprocal connection probability measured across all trials:
    print("  measured reciprocal probability: %0.4f±%0.4f" % (np.mean(rprobs), np.std(rprobs)))
    # Null-hypothesis estimated reciprocal connection probability measured across all trials (Pc^2)
    # (this is the value that the experimenter would compare to their measured reciprocal connection probability) 
    print("  estimatd reciprocal probability: %0.4f±%0.4f" % (np.mean(ex_rprobs), np.std(ex_rprobs)))
 
    # Ratio of measured versus estimated reciprocal connection probability. A value > 1 is
    # meant to indicate that reciprocal connections are "preferred" above random chance
    print("  measured/estimated ratio: %0.4f±%0.4f" % (np.mean(ratios), np.std(ratios)))

    return results


if __name__ == '__main__':
    import pyqtgraph as pg

    pg.mkQApp()

    params = pg.parametertree.Parameter.create(name='params', type='group', children=[
        {'name': 'n_types', 'type': 'int', 'value': 2},
        {'name': 'n_cells', 'type': 'int', 'value': 6},
        {'name': 'n_expts', 'type': 'int', 'value': 10},
        {'name': 'n_trials', 'type': 'int', 'value': 200},
        {'name': 'run', 'type': 'action'},
    ])

    vs = pg.QtWidgets.QSplitter(pg.QtCore.Qt.Vertical)

    pt = pg.parametertree.ParameterTree()
    vs.addWidget(pt)
    pt.setParameters(params)

    wab_table = pg.QtWidgets.QTableWidget()
    vs.addWidget(wab_table)

    hs = pg.QtWidgets.QSplitter(pg.QtCore.Qt.Horizontal)
    hs.addWidget(vs)

    gl = pg.GraphicsLayoutWidget()
    hs.addWidget(gl)

    c_plt = gl.addPlot(0, 0, labels={'bottom': 'measured connection probability'})
    r_plt = gl.addPlot(0, 1, labels={'bottom': 'measured reciprocal probability'})
    rr_plt = gl.addPlot(1, 0, labels={'bottom': 'cp^2', 'left':'measured reciprocal probability'})
    # rr_plt.setAspectLocked(1)
    ratio_plt = gl.addPlot(1, 1, labels={'bottom': 'measured reciprocal / cp^2'})
    hs.show()

    def run():
        global results  # allow inspection from console
        n_types = params['n_types']
        Wab = np.array([[float(str(wab_table.item(i, j).text())) for j in range(n_types)] for i in range(n_types)])
        with pg.BusyCursor():
            results = run_expt(Wab, n_cells=params['n_cells'], n_expts=params['n_expts'], n_trials=params['n_trials'])

        cprobs = results['conn'].sum(axis=1) / results['probed'].sum(axis=1)
        rprobs = results['recip'].sum(axis=1) / results['probed'].sum(axis=1)
        ex_rprobs = cprobs**2
        ratios = rprobs / ex_rprobs
        
        y = pg.pseudoScatter(cprobs)
        c_plt.plot(cprobs, y, clear=True, pen=None, symbol='o')
        c_plt.addLine(x=Wab.mean())
        

        y = pg.pseudoScatter(rprobs)
        r_plt.plot(rprobs, y, clear=True, pen=None, symbol='o')

        y = pg.pseudoScatter(ratios)
        ratio_plt.plot(ratios, y, clear=True, pen=None, symbol='o')
        ratio_plt.addLine(x=1)
        
        rr_plt.plot(ex_rprobs, rprobs, clear=True, pen=None, symbol='o')
        l = pg.InfiniteLine(angle=45)
        rr_plt.addItem(l)

    def set_rc():
        n_types = params['n_types']
        wab_table.setRowCount(n_types)
        wab_table.setColumnCount(n_types)
        for i in range(n_types):
            for j in range(n_types):
                item = pg.QtWidgets.QTableWidgetItem('0.1')
                wab_table.setItem(i, j, item)
            wab_table.setColumnWidth(i, 40)
            wab_table.setRowHeight(i, 40)
        
    set_rc()
    run()
    
    params.child('run').sigActivated.connect(run)
    params.child('n_types').sigValueChanged.connect(set_rc)
    


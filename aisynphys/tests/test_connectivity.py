#%%
import pytest
import numpy as np
import pdb
from aisynphys.cell_class import CellClass
from aisynphys.connectivity import CorrectionModel, GaussianModel, BinaryModel, ErfModel, ei_correct_connectivity



def data_snippets(conntype='finite'):
    # small set of pairs from the data (values are approximated)
    keys = ['pre_axon_length', 'avg_pair_depth', 'detection_power', 'lateral_distance', 'has_synapse', 'pre_cell', 'n_in_test_spikes']
    metrics = {}
    metrics[keys[0]] = np.array([2e-4, 1e-4, 2e-4, None, None])
    metrics[keys[1]] = np.array([3.3e-5, 9.0e-5, 7.7e-5, 1.27e-4, 1.1e-4])
    metrics[keys[2]] = np.array([4.57, 4.52, np.nan, 4.27, 4.46])
    metrics[keys[3]] = np.array([5.3e-5, 2.7e-5, 1.99e-4, 4.3e-5, 6.7e-5])
    if conntype == 'finite':
        metrics[keys[4]] = np.array([False, False, False, False, True])
    elif conntype == 'zero':
        metrics[keys[4]] = np.array([False]*5)
    elif conntype == 'full':
        metrics[keys[4]] = np.array([True]*5)
        
    # pre_cell type
    cell = lambda: None # initiating a function so that I can add an attribute to it
    cell.cell_class_nonsynaptic = 'in'
    metrics[keys[5]] = [cell] * 5
    metrics[keys[6]] = [800] * 5

    pairs = []
    for i in range(5):
        pair = lambda: None
        for k in keys:
            setattr(pair, k, metrics[k][i])
        pairs.append(pair)
    
    variables = [metrics[keys[i]] for i in [3, 0, 1, 2]]
    conn = metrics[keys[4]]

    return pairs, variables, conn

def test_correction_model():
    # test preparations
    ei_classes = {'excit': CellClass(cell_class_nonsynaptic='ex', name='ex'), 'inhib': CellClass(cell_class_nonsynaptic='in', name='in')}
    correction_metrics = {
        'lateral_distance': {'model': GaussianModel, 'init': (0.1, 100e-6), 'bounds': ((0.001, 1), (100e-6, 100e-6))}, 
        'pre_axon_length': {'model': BinaryModel, 'init': (0.1, 200e-6, 0.5), 'bounds': ((0.001, 1), (200e-6, 200e-6), (0.1, 0.9))}, 
        'avg_pair_depth': {'model': ErfModel, 'init': (0.1, 30e-6, 30e-6), 'bounds': ((0.01, 1), (10e-6, 200e-6), (-100e-6, 100e-6))},
        'detection_power': {'init': (0.1, 1.0, 3.0), 'model': ErfModel, 'bounds': ((0.001, 1), (0.1, 5), (2, 5)),
                              'constraint': (0.6745, 4.6613)}, 
    }
    
    # test finite connection cases
    pairs, variables, conn = data_snippets()
    corr_model = ei_correct_connectivity(ei_classes, correction_metrics, pairs)
    # fixing the correction parameters from a larger pool of data
    correction_parameters = [[np.array([0.1, 100e-6]),                   # lateral_distance (Gaussian)
                              np.array([6.939e-2, 2.000e-4, 6.003e-1]),  # pre_axon_length (Binary)
                              np.array([6.893e-2, 4.021e-5, 2.848e-5]),  # avg_pair_depth (Erf)
                              np.array([9.847e-2, 6.845e-1, 4.200])],    # detection_power (Erf)
                             [np.array([0.1, 100e-6]),
                              np.array([1.227e-1, 2.000e-4, 3.634e-1]),
                              np.array([1.095e-1, 4.130e-5, 2.700e-5]),
                              np.array([1.967e-1, 4.396e-1, 4.365])]]
    
    corr_model.correction_parameters = correction_parameters
    
    assert corr_model.pmax == 0.1

    #fit = CorrectionModel.fit(corr_model, variables, conn, excinh=1)
    fit = corr_model.fit(variables, conn, excinh=1)
    assert corr_model.pmax == pytest.approx(0.74590, 1e-4)
    assert fit.x == pytest.approx(0.74590, 1e-4)
    assert corr_model.likelihood(variables, conn) == pytest.approx(-1.95075, 1e-4)
    assert fit.fun == pytest.approx(1.95075, 1e-4)
    
    assert fit.cp_ci[0] == corr_model.pmax
    assert fit.cp_ci[1] == pytest.approx(0.0478545, 1e-4) # lower MINOS value
    assert fit.cp_ci[2] == pytest.approx(2.17707, 1e-4) # upper MINOS value

    # test zero-connection cases
    pairs, variables, conn = data_snippets('zero')
    fit = CorrectionModel.fit(corr_model, variables, conn, excinh=1)
    assert fit.cp_ci[0] == pytest.approx(0, abs=1e-6)
    assert fit.cp_ci[1] == pytest.approx(0, abs=1e-6) # lower MINOS value
    assert fit.cp_ci[2] == pytest.approx(1.03155, 1e-4) # upper MINOS value

    # test full-cnnection cases
    pairs, variables, conn = data_snippets('full')
    fit = CorrectionModel.fit(corr_model, variables, conn, excinh=1)
    assert fit.cp_ci[0] == pytest.approx(3, 1e-4) # maxed out
    assert fit.cp_ci[1] == pytest.approx(1.87319, 1e-4) # lower MINOS value
    assert fit.cp_ci[2] == pytest.approx(3, 1e-4) # upper MINOS value

    return
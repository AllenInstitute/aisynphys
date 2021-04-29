# -*- coding: utf-8 -*- 
import warnings
from collections import OrderedDict
import numpy as np
from statsmodels.stats.proportion import proportion_confint
import scipy.optimize
from scipy.special import erf
from .util import optional_import
from aisynphys.database import default_db as db
iminuit = optional_import('iminuit')


def connectivity_profile(connected, distance, bin_edges):
    """
    Compute connection probability vs distance with confidence intervals.

    Parameters
    ----------
    connected : boolean array
        Whether a synaptic connection was found for each probe
    distance : array
        Distance between cells for each probe
    bin_edges : array
        The distance values between which connections will be binned

    Returns
    -------
    xvals : array
        bin edges of returned connectivity values
    prop : array
        connected proportion in each bin
    lower : array
        lower proportion confidence interval for each bin
    upper : array
        upper proportion confidence interval for each bin

    """
    mask = np.isfinite(connected) & np.isfinite(distance)
    connected = connected[mask]
    distance = distance[mask]

    n_bins = len(bin_edges) - 1
    upper = np.zeros(n_bins)
    lower = np.zeros(n_bins)
    prop = np.zeros(n_bins)
    for i in range(n_bins):
        minx = bin_edges[i]
        maxx = bin_edges[i+1]

        # select points inside this window
        mask = (distance >= minx) & (distance < maxx)
        pts_in_window = connected[mask]
        # compute stats for window
        n_probed = pts_in_window.shape[0]
        n_conn = pts_in_window.sum()
        if n_probed == 0:
            prop[i] = np.nan
            lower[i] = 0
            upper[i] = 1
        else:
            prop[i] = n_conn / n_probed
            ci = connection_probability_ci(n_conn, n_probed)
            lower[i] = ci[0]
            upper[i] = ci[1]

    return bin_edges, prop, lower, upper


def measure_distance(pair_groups, window):
    """Given a description of cell pairs grouped together by cell class,
    return a structure that describes connectivity as a function of distance between cell classes.
    
    Parameters
    ----------
    pair_groups : OrderedDict
        Output of `cell_class.classify_pairs`
    window: float
        binning window for distance
    """

    results = OrderedDict()
    for key, class_pairs in pair_groups.items():
        pre_class, post_class = key

        connected, distance = pair_distance(class_pairs, pre_class) 
        bin_edges = np.arange(0, 500e-6, window)
        xvals, cp, lower, upper = connectivity_profile(connected, distance, bin_edges)

        results[(pre_class, post_class)] = {
        'bin_edges': bin_edges,
        'conn_prob': cp,
        'lower_ci': lower,
        'upper_ci': upper,
        }

    return results


def pair_distance(class_pairs, pre_class):
    """Given a list of cell pairs return an array of connectivity and distance for each pair.
    """

    connected = []
    distance = []

    for pair in class_pairs:
        probed = pair_was_probed(pair, pre_class.output_synapse_type)
        if probed and pair.distance is not None:
            connected.append(pair.has_synapse)
            distance.append(pair.distance)

    connected = np.asarray(connected).astype(float)
    distance = np.asarray(distance)

    return connected, distance


def measure_connectivity(pair_groups, alpha=0.05, sigma=None, fit_model=None, correction_model=None, dist_measure='distance'):
    """Given a description of cell pairs grouped together by cell class,
    return a structure that describes connectivity between cell classes.
    
    Parameters
    ----------
    pair_groups : OrderedDict
        Output of `cell_class.classify_pairs`
    alpha : float
        Alpha value setting confidence interval width (default is 0.05)
    sigma : float | None
        Sigma value for distance-adjusted connectivity (see 
        ``distance_adjysted_connectivity()``). If None, then adjusted
        values are omitted from the result.
    fit_model : ConnectivityModel | None
        ConnectivityModel subclass to fit Cp vs distance profile. If combined with
        sigma the fit will be fixed to that sigma. If None, then fit results
        are ommitted from the results
    correction_model: ConnectivityModel | None
        ConnectivityModel subclass used to correct fit_model based on other metrics 
        which are set within ConnectivityModel.variables
    dist_measure : str
        Which distance measure to use when calculating connection probability.
        Must be one of 'distance', 'lateral_distance', 'vertical_distance' columns
        from Pair table in SynPhys database

    Returns
    -------
    result : dict
        Keys are the same as in the *pair_groups* argument. Values are dictionaries
        containing connectivity results for each pair group::
        
            {n_probed=int, n_connected=int, probed_pairs=list, connected_pairs=list,
             connection_probability=(cp, lower_ci, upper_ci), 
             adjusted_connectivity=(cp, lower_ci, upper_ci)}
    """    
    results = OrderedDict()
    for key, class_pairs in pair_groups.items():
        pre_class, post_class = key
        
        probed_pairs = [p for p in class_pairs if pair_was_probed(p, pre_class.output_synapse_type)]
        connections_found = [p for p in probed_pairs if p.has_synapse]
        
        gaps_probed = [p for p in class_pairs if pair_probed_gj(p)]
        gaps_found = [p for p in gaps_probed if p.has_electrical]

        n_connected = len(connections_found)
        n_probed = len(probed_pairs)
        n_gaps_probed = len(gaps_probed)
        n_gaps = len(gaps_found)
        conf_interval_cp = connection_probability_ci(n_connected, n_probed, alpha=alpha)
        conn_prob = float('nan') if n_probed == 0 else n_connected / n_probed
        conf_interval_gap = connection_probability_ci(n_gaps, n_gaps_probed, alpha=alpha)
        gap_prob = float('nan') if n_gaps_probed == 0 else n_gaps / n_gaps_probed

        results[(pre_class, post_class)] = {
            'n_probed': n_probed,
            'n_connected': n_connected,
            'n_gaps_probed': n_gaps_probed,
            'n_gaps': n_gaps,
            'connection_probability': (conn_prob,) + conf_interval_cp,
            'gap_probability': (gap_prob,) + conf_interval_gap,
            'connected_pairs': connections_found,
            'gap_pairs': gaps_found,
            'probed_pairs': probed_pairs,
        }

        if sigma is not None or fit_model is not None:
            distances = np.array([getattr(p, dist_measure) for p in probed_pairs], dtype=float)
            connections = np.array([p.synapse for p in probed_pairs], dtype=bool)
            gap_distances = np.array([getattr(p, dist_measure) for p in gaps_probed], dtype=float)
            gaps = np.array([p.has_electrical for p in gaps_probed], dtype=bool)
            mask = np.isfinite(distances) & np.isfinite(connections)
            results[(pre_class, post_class)]['probed_distances'] = distances[mask]
            results[(pre_class, post_class)]['connected_distances'] = connections[mask]
            mask2 = np.isfinite(gap_distances) & np.isfinite(gaps)
            results[(pre_class, post_class)]['gap_probed_distances'] = gap_distances[mask2]
            results[(pre_class, post_class)]['gap_distances'] = gaps[mask2]
            fit = fit_model.fit(distances[mask], connections[mask], method='L-BFGS-B', fixed_size=sigma)
            results[(pre_class, post_class)]['connectivity_fit'] = fit
            gap_fit = fit_model.fit(gap_distances[mask2], gaps[mask2], method='L-BFGS-B', fixed_size=sigma)
            results[(pre_class, post_class)]['gap_fit'] = gap_fit
        if sigma is not None:
            adj_conn_prob, adj_lower_ci, adj_upper_ci = distance_adjusted_connectivity(distances[mask], connections[mask], sigma=sigma, alpha=alpha)
            results[(pre_class, post_class)]['adjusted_connectivity'] = (adj_conn_prob, adj_lower_ci, adj_upper_ci)
            adj_gap_junc, adj_lower_gj_ci, adj_upper_gj_ci = distance_adjusted_connectivity(gap_distances[mask2], gaps[mask2], sigma=sigma, alpha=alpha)
            results[(pre_class, post_class)]['adjusted_gap_junction'] = (adj_gap_junc, adj_lower_gj_ci, adj_upper_gj_ci)
        if correction_model is not None:
            # Here it performs corrected p_max fit if there are relevant variables.
            if hasattr(correction_model, 'correction_variables'): # correction model.
                correction_model.size = sigma # override it
                variables = [distances[mask]]
                for variable in correction_model.correction_variables:
                    var_extract = np.array([getattr(p, variable) for p in probed_pairs], dtype=float)
                    variables.append(var_extract[mask])
                if len(class_pairs) == 0: # empty list
                    excinh = 0 # doesn't matter which, so give exc.
                elif class_pairs[0].pre_cell.cell_class_nonsynaptic == 'ex':
                    excinh = 0
                else:
                    excinh = 1

                corr_fit = correction_model.fit(correction_model, variables, connections[mask], excinh=excinh)
                if mask.sum() == 0: # no probing
                    corr_fit.x = np.nan # needed not to report the initial value
                results[(pre_class, post_class)]['connectivity_correction_fit'] = corr_fit
                # for gap junctions, this analysis won't be relevant, so I won't assign gap_fit for now.              
    
    return results


def connection_probability_ci(n_connected, n_probed, alpha=0.05):
    """Return confidence intervals on the probability of connectivity, given the
    number of putative connections probed vs the number of connections found.
    
    Currently this simply calls `statsmodels.stats.proportion.proportion_confint`
    using the "beta" method.
    
    Parameters
    ----------
    n_connected : int
        The number of observed connections in a sample
    n_probed : int
        The number of probed (putative) connections in a sample; must be >= n_connected
        
    Returns
    -------
    lower : float
        The lower confidence interval
    upper : float
        The upper confidence interval
    """
    assert n_connected <= n_probed, "n_connected must be <= n_probed"
    if n_probed == 0:
        return (0, 1)
    return proportion_confint(n_connected, n_probed, alpha=alpha, method='beta')


def pair_was_probed(pair, synapse_type):
    """Return boolean indicating whether a cell pair was "probed" for either 
    excitatory or inhibitory connectivity.
    
    Currently this is determined by checking that either `n_ex_test_spikes` or 
    `n_in_test_spikes` is greater than 10, depending on the value of *synapse_type*.
    This is an arbitrary limit that trades off between retaining more data and
    rejecting experiments that were not adequately sampled. 
    
    Parameters
    ----------
    synapse_type : str
        Must be either 'ex', 'in', or None. If None, then the pair is considered probed
        if it passes criteria for _both_ 'ex' and 'in'.
    """
    test_spike_limit = 10

    assert synapse_type in ('ex', 'in', None), "synapse_type must be 'ex', 'in', or None"
    if synapse_type is None:
        return (pair.n_ex_test_spikes > test_spike_limit) and (pair.n_in_test_spikes > test_spike_limit)
    else:
        qc_field = 'n_%s_test_spikes' % synapse_type
        return getattr(pair, qc_field) > test_spike_limit


def pair_probed_gj(pair):
    """Return boolean indicateing whether a cell pair was "probed" for a gap junction.

    Checks both the "pre" and "post" synaptic cells for long-pulse stimuli from which
    gap junctions were identified and quantified.
    """
    
    pre_electrode = pair.pre_cell.electrode
    post_electrode = pair.post_cell.electrode
    pre_stims = set([rec.stim_name for rec in pre_electrode.recordings])
    post_stims = set([rec.stim_name for rec in post_electrode.recordings])
    return any('TargetV' in s for s in pre_stims) and any('TargetV' in s for s in post_stims)


def distance_adjusted_connectivity(x_probed, connected, sigma, alpha=0.05):
    """Return connectivity and binomial confidence interval corrected for the distances
    at which connections are probed.
    
    This function models connectivity as a gaussian curve with respect to intersomatic 
    distance; two cells are less likely to be connected if the distance between them is
    large. Due to this relationship between distance and connectivity, simple measures
    of connection probability are sensitive to the distances at which connections are
    tested. This function returns a connection probability and CI that are adjusted
    to normalize for these distances.
    
    The returned *pmax* is also the maximum value of a gaussian connectivity
    profile that most closely matches the input data using a maximum likelihood
    estimation. In an ideal scenario we would use this to determine both
    the gaussian _amplitude_ and _sigma_ that most closely match the data. In practice,
    however, very large N is required to constrain sigma so we instead use a fixed sigma
    and focus on estimating the amplitude. Keep in mind that the corrections performed
    here rely on the assumption that connectivity falls off with distance as a gaussian
    of this fixed sigma value; use caution interpreting these results in cases where
    that assumption may be wrong (for example, for inter-laminar connectivity in cortex).

    Another practical constraint is that exact estimates and confidence intervals are 
    expensive to compute, so we use a close approximation that simply scales the 
    connectivity proportion and binomial confidence intervals using the average connection
    probability from the gaussian profile. 

    For more detail on how this method was chosen and developed, see
    aisynphys/doc/connectivity_vs_distance.ipynb.

    Parameters
    ----------
    x_probed : array
        Array containing intersomatic distances of pairs that were probed for connectivity.
    connected : bool array
        Boolean array indicating which pairs were connected.
    sigma : float
        Gaussian sigma value defining the width of the connectivity profile to fit to *x_probed* and *connected*.
    alpha : float
        Alpha value setting the width of the confidence interval. Default is 0.05, giving a 95% CI.

    Returns
    -------
    pmax : float
        Maximum probability value in the gaussian profile (at x=0)
    lower : float
        Lower edge of confidence interval
    upper : float
        Upper edge of confidence interval

    """
    # mean connection probability of a gaussian profile where cp is 1.0 at the center,
    # sampled at locations probed for connectivity
    mean_cp = np.exp(-x_probed**2 / (2 * sigma**2)).mean()

    n_conn = connected.sum()
    n_test = len(x_probed)

    # estimated pmax is just the proportion of connections multiplied by a scale factor
    est_pmax = (n_conn / n_test) / mean_cp
    
    # and for the CI, we can just use a standard binomial confidence interval scaled by the same factor
    lower, upper = connection_probability_ci(n_conn, n_test, alpha)
    
    return est_pmax, lower / mean_cp, upper / mean_cp


class ConnectivityModel:
    """Base class for modeling and fitting connectivity vs intersomatic distance data.

    Implements fitting by maximum likelihood estimation.
    """
    def generate(self, x, seed=None):
        """Generate a random sample of connectivity, given distances at which connections are tested.
        
        The *seed* parameter allows a random seed to be provided so that results are frozen across executions.
        """
        p = self.connection_probability(x)
        rng = np.random.RandomState(seed)
        return rng.random(size=len(p)) < p

    def likelihood(self, x, conn):
        """Log-likelihood for maximum likelihood estimation

        LLF = Σᵢ(𝑦ᵢ log(𝑝(𝐱ᵢ)) + (1 − 𝑦ᵢ) log(1 − 𝑝(𝐱ᵢ)))
        """
        assert np.issubdtype(conn.dtype, np.dtype(bool))
        p = self.connection_probability(x)
        return np.log(p[conn]).sum() + np.log((1-p)[~conn]).sum()

    def connection_probability(self, x):
        """Return the probability of seeing a connection from cell A to cell B given 
        the intersomatic distance *x*.
        
        This does _not_ include the probability of the reverse connection B→A.
        """
        raise NotImplementedError('connection_probability must be implemented in a subclass')
    
    @classmethod
    def err_fn(cls, params, *args):
        model = cls(*params)
        return -model.likelihood(*args)

    @classmethod
    def fit(cls, x, conn, init=(0.1, 150e-6), bounds=((0.001, 1), (10e-6, 1e-3)), fixed_size=None, fixed_max=None, **kwds):
        n = 6
        p_bins = np.linspace(bounds[0][0], bounds[0][1], n)
        if fixed_max is not None:
            p_bins = [fixed_max]*2
        s_bins = np.linspace(bounds[1][0], bounds[1][1], n)
        if fixed_size is not None:
            s_bins = [fixed_size]*2
        best = None
        # Most minimization methods fail to find the global minimum for this problem.
        # Instead, we systematically search over a large (n x n) range of the parameter space
        # and pick the best overall fit.
        for p_bin in zip(p_bins[:-1], p_bins[1:]):
            for s_bin in zip(s_bins[:-1], s_bins[1:]):
                init = (np.mean(p_bin), np.mean(s_bin))
                bounds = (p_bin, s_bin)

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")    
                    fit = scipy.optimize.minimize(
                        cls.err_fn, 
                        x0=init, 
                        args=(x, conn),
                        bounds=bounds,
                        **kwds,
                    )
                    if best is None or fit.fun < best.fun:
                        best = fit
        
        ret = cls(*best.x)
        ret.fit_result = best
        return ret


class SphereIntersectionModel(ConnectivityModel):
    """Model connection probability as proportional to the volume overlap of intersecting spheres.

    Parameters
    ----------
    pmax : float
        Maximum connection probability (at 0 intersomatic distance)
    size : float
        Radius of a single cell
    density : float
        Average number of synapses per m^3 connecting one cell to another, within their overlapping volume.        

    Note: the *pmax* and *density* parameters are redundant; can be instantiated using either 
    (pmax, size) or (density, size).
    """
    def __init__(self, pmax=None, size=None, density=None):
        assert size is not None
        assert (pmax is not None) or (density is not None)
            
        self._pmax = pmax
        self.size = size
        self._density = density

    @property
    def r(self):
        """Radius of cell; overlapping region extends to 2*r.
        """
        return self.size
        
    @property
    def density(self):
        """Average density of synapses in overlapping regions (synapses / m^3)"""
        if self._density is None:
            # calculate synapse density needed to get pmax
            return np.log(1 - self.pmax) / (-(4/3) * np.pi * self.r**3)
        else:
            return self._density

    @density.setter
    def density(self, d):
        self._density = d
        self._pmax = None

    @property
    def pmax(self):
        """Maximum connection probability at intersomatic distance=0
        """
        if self._pmax is None:
            return self.connection_probability(0)
        else:
            return self._pmax

    @pmax.setter
    def pmax(self, p):
        self._pmax = p
        self._density = None
    
    def connection_probability(self, x):
        # In a Poisson process, the probability of seeing 0 events in an interval is exp(-λ), where
        # λ is the average number of events in an interval. If we let *density* in this case be the average
        # number of synapses per unit volume, then λ = density * V, and the probability of seeing no events
        # in V is just P(V,density) = exp(-denstiy*V).
        
        r = self.r
        density = self.density
        
        # volume of overlap between spheres at distance x
        v = self.volume_overlap(x)
        
        # probability of 1 or more synapse given v and density
        p = 1 - np.exp(-v * density)
        
        # clip probability at 0.005 to allow for a small number of
        # connections past 2*r (otherwise fitting becomes impossible)
        return np.clip(p, 0.005, 1)
    
    def volume_overlap(self, x):
        r = self.r
        v = (4 * np.pi / 3) * r**3 - np.pi * x * (r**2 - x**2/12)
        return np.where(x < r*2, v, 0)


class ExpModel(ConnectivityModel):
    """Model connection probability as an exponential decay
    
    Parameters
    ----------
    pmax : float
        Maximum connection probability (at 0 intersomatic distance)
    size : float
        Exponential decay tau
    
    """
    def __init__(self, pmax, size):
        self.pmax = pmax
        self.size = size
    
    def connection_probability(self, x):
        return self.pmax * np.exp(-x / self.size)

    
class LinearModel(ConnectivityModel):
    """Model connection probability as linear decay to 0

    Parameters
    ----------
    pmax : float
        Maximum connection probability (at 0 intersomatic distance)
    size : float
        Distance at which connection probability reaches its minimum value
        
    Note: in this model, connection probability goes only to 0.5% at minimum.
    If the minimum is allowed to go to 0, then fitting becomes more difficult.
    """
    def __init__(self, pmax, size):
        self.pmax = pmax
        self.size = size

    def connection_probability(self, x):
        return np.clip(self.pmax * (self.size - x) / self.size, 0.005, 1)


class GaussianModel(ConnectivityModel):
    """Model connection probability as gaussian centered at x=0

    Parameters
    ----------
    pmax : float
        Maximum connection probability (at 0 intersomatic distance)
    size : float
        Gaussian sigma
    """
    def __init__(self, pmax, size):
        self.pmax = pmax
        self.size = size

    def connection_probability(self, x):
        return self.pmax * np.exp(-x**2 / (2 * self.size**2))

    @staticmethod
    def correction_func(params, x):
        return np.exp(-x**2 / (2 * params[1]**2))

class CorrectionModel(ConnectivityModel):
    """ Connectivity model with corrections for potential biases

    Parameters
    ----------
    pmax : float
        Maximum connection probability (at 0 intersomatic distance)
    correction_variables : list of strings
        Names of correction variables.
        These names should be defined in each Pair in pair_groups when measure_connectivity is called.
    correction_functions : list of functions
        Functions used to correct estimating connection probability.
        The following formats need to be valid to execute properly.
        correction_functions[i](correction_parameters[excinh][i], pair[correction_variables[i]])
        where is a Pair instance in pair_group (1st argument of measure_connectivity)
    correction_parameters : list of list of [array or [list of float]]
        Parameters used aside with correction_variables. Fixed during the fit.
        The first index is 0 (pre-synaptic excitatory) or 1 (pre-synaptic inhibitory)
    do_minos : bool
        Use MINOS algorithm for estimating the 95% confidence interval (default: True).
        If False, it uses the same approximation method used for the distance adjustment.
        
    Note: correction_variables, correction_parameters, and correction_functions should have the same lengths.
    """
    def __init__(self, pmax, correction_variables, correction_functions, correction_parameters, do_minos=True):
        self.pmax = pmax
        self.correction_variables = correction_variables # list of strings (names of correction variables)
        self.correction_functions = correction_functions # list of functions
        self.correction_parameters = correction_parameters # list of list of parameters for correction functions
        self.do_minos = do_minos


    def connection_probability(self, x):
        # x is expected to be [distance, correction_var1, correction_var2,...].
        # the final probability is modeled as a product of multiple corrections.
        correction = 1.0
        for i in range(len(x)):
            # replace None with nan to make the following code to work
            v = np.array([np.nan if el is None else el for el in x[i]])
            corrval = self.correction_functions[i](self.correction_parameters[self.excinh][i], v)
            correction *= np.nan_to_num(corrval, nan=1.0)
        return np.clip(self.pmax * correction, 0.0, 1.0)


    def nll(self, pmax, x, conn):
        self.pmax = pmax # override existing value
        return -self.likelihood(x, conn)


    def fit(self, x, conn, init=(0.1), bounds=((0.0, 3.0)), excinh=None, **kwds):
        self.excinh = excinh # setting the cell class...
        fit = iminuit.minimize(
                        self.nll,
                        x0=init, 
                        args=(x, conn),
                        bounds=bounds,
                        **kwds,
                    )

        if len(conn) == 0:
            fit.cp_ci = (np.nan, np.nan, np.nan)
            return fit

        cp = fit.x
        import pdb
        
        if self.do_minos: # MINOS (Likelihood-based CI estimation)
            # do minuit calculation
            cl = 0.95

            try:
                fit.minuit.minos(cl=cl)
            except RuntimeError as e: # Catch RuntimeError from iminuit
                if conn.sum() == 0:
                    # zero connection case. special treatment
                    fit.x = 0
                    lower, upper = manual_ci_search(fit, cl, zero_connection=True)
                    fit.cp_ci = (cp, lower, upper)
                    pdb.set_trace()
                    return fit
                else:
                    lower, upper = manual_ci_search(fit, cl)
                    fit.cp_ci = (cp, lower, upper)
                    print(e)
                    return fit
            except Exception as e:
                # re-raise if unexpected exception is observed
                print(e)
                raise e

            cp = fit.minuit.params['x0'].value
            lower = cp  + fit.minuit.merrors['x0'].lower
            upper = cp  + fit.minuit.merrors['x0'].upper
            if fit.minuit.merrors['x0'].is_valid == False:
                # do it once more (sometimes helps)
                fit.minuit.minos(cl=cl)
                if fit.minuit.merrors['x0'].is_valid == False:
                    # it it doesn't work, manually search the value
                    lower, upper = manual_ci_search(fit, cl)

        else: # Approximation method
            self.pmax = 1.0
            mean_adjustment = self.connection_probability(x).mean()
            # to calculate CI, get the adjustment values from the instance.
            n_conn = conn.sum()
            n_test = len(conn)
            est_pmax = (n_conn / n_test) / mean_adjustment
        
            # and for the CI, we can just use a standard binomial confidence interval scaled by the same factor
            lower, upper = connection_probability_ci(n_conn, n_test)

            lower = lower / mean_adjustment
            upper = upper / mean_adjustment
            fit.cp_ci = (est_pmax, lower, upper)

        fit.cp_ci = (cp, lower, upper)

        return fit


def manual_ci_search(fit, cl, zero_connection=False):
    """Return confidence intervals on the probability of connectivity, given the
    failed execution of the minuit fit.
    
    Currently this simply calls `statsmodels.stats.proportion.proportion_confint`
    using the "beta" method.
    
    Parameters
    ----------
    fit : Minuit
        Result of the fit returned by iminuit
    cl : float
        Confidence level for the CI estimate
    zero_connection : bool
        If true, set return 0 as the lower CI bound, and only evaluate the upper CI bound
        
    Returns
    -------
    lower : float
        The lower confidence interval
    upper : float
        The upper confidence interval
    """
    from scipy.stats import chi2
    factor = chi2(1).ppf(cl)
    fitfun = lambda x: (fit.minuit.fcn([x]) - fit.minuit.fcn([fit.x]) - factor * fit.minuit.errordef)**2
    if zero_connection:
        lower = 0.0
        upper = iminuit.minimize(fitfun, 0.001, bounds=[[fit.x, 5]])
        return (lower, upper.x)
    else:
        lower = iminuit.minimize(fitfun, fit.x * 0.9, bounds=[[0, fit.x]])
        upper = iminuit.minimize(fitfun, fit.x * 1.1, bounds=[[fit.x, 5]])
        return (lower.x, upper.x)
    


class ErfModel(ConnectivityModel):
    """Model connection probability correction as error function

    Parameters
    ----------
    pmax : float
        Maximum connection probability (at 0 intersomatic distance)
    size : float
        Sigma of the integrand Gaussian of the error function
    midpoint : float
        Mu (center) of the integrand Gaussian of the error function

    Note: You can give a constraint when this model is fit.
    Pass a constrant argument (2-element tuple with (sigma_multiplier, stop_point)) to fit function.
    The actual constraint equation is the following format.
        midpoints + constraint[0] * size < constraint[1]
    To make quartile saturation at quartile of the data point of detection power, use the following parameters.
    constraint = (0.6745, 4.5655)
    See Connectivity Corrections Fit.ipynb for more details.
    """
    def __init__(self, pmax, size, midpoint):
        self.pmax = pmax
        self.size = size
        self.midpoint = midpoint
    
    def connection_probability(self, x):
        return self.pmax * 0.5 * (1.0 + erf((x - self.midpoint) / (np.sqrt(2) * self.size)))

    @staticmethod
    def correction_func(params, x):
        return 0.5 * (1.0 + erf((x - params[2]) / (np.sqrt(2) * params[1])))

    @classmethod
    def fit(cls, x, conn, init, bounds, constraint=None, fixed_size=False, fixed_max=False):
        # constraint should be a 2-element tuple with (sigma_multiplier, stop_point)
        if constraint == None:
            fit = iminuit.minimize(
                cls.err_fn,
                x0=init,
                args=(x, conn),
                bounds=bounds,
            )
        else:
            # midpoints + constraint[0] * size < constraint[1]
            # to make it quartile of the detection power, use the following parameters
            # constraint[0] == 0.6745 (specify quartile point using sigma)
            # constraint[1] == 4.5655 observed value in the data (may change if the data change)
            con_array = np.array([[0, constraint[0], 1]]) # 1 x 3 matrix
            constraints = scipy.optimize.LinearConstraint(con_array, lb=-np.inf, ub=constraint[1])
            fit = scipy.optimize.minimize(
                cls.err_fn,
                x0=init,
                args=(x, conn),
                bounds=bounds,
                constraints=constraints
            )
        ret = cls(*fit.x)
        ret.fit_result = fit
        return ret


class BinaryModel(ConnectivityModel):
    """Model connection probability correction as binary function

    Parameters
    ----------
    pmax : float
        Maximum connection probability (at 0 intersomatic distance)
    size : float
        Threshold for applying adjustment. If variable < size, adjustment is applied.
    adjustment : float
        Ratio between adjusted and non-adjusted probability.
        Resulting probability after adjustment is pmax * adjustment
    """
    def __init__(self, pmax, size, adjustment):
        self.pmax = pmax
        self.size = size
        self.adjustment = adjustment
    
    def connection_probability(self, x):
        val = np.ones_like(x)
        val[x < self.size] *= self.adjustment
        return self.pmax * val

    @staticmethod
    def correction_func(params, x):
        val = np.ones_like(x)
        val[x < params[1]] *= params[2]
        return val

    @classmethod
    def fit(cls, x, conn, init, bounds, fixed_size=False, fixed_max=False):
        fit = iminuit.minimize(
            cls.err_fn,
            x0=init,
            args=(x, conn),
            bounds=bounds,
        )
        ret = cls(*fit.x)
        ret.fit_result = fit
        return ret


class FixedSizeModelMixin:
    @classmethod
    def fit(cls, x, conn, init=0.1, bounds=(0.001, 1), **kwds):
        fixed_size = kwds.pop('size', 130e-6)
        
        if conn.sum() == 0:
            # fitting falls apart at 0; just return the obvious result
            return cls(0, fixed_size)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    
            fit = scipy.optimize.minimize(
                cls.err_fn, 
                x0=(init,), 
                args=(fixed_size, x, conn),
                bounds=(bounds,),
                **kwds,
            )
        
        ret = cls(fit.x[0], fixed_size)
        ret.fit_result = fit
        return ret
    
    @classmethod
    def err_fn(cls, params, *args):
        (pmax,) = params
        size, x, conn = args
        model = cls(pmax, size)
        return -model.likelihood(x, conn)


def recip_connectivity_profile(probes_1, probes_2, bin_edges):
    """
    Given probed pairs of two pair groups what is the normalized probability of finding a reciprocal cnx 
    
    Parameters
    ----------
    probes_1 : list of Pair
        probed pairs for first connection type (A->B)
    probes_2 : list of Pair
        probed pairs for second connection type (B->A)
    bin_edges : array
        The distance values between which connections will be binned

    Output
    ------
    norm_cp_r : array
        distance binned reciprocal connection probability normalized to cp_1 * cp_2
    recip_conn : boolean array
        whether or not a reciprocal connection exists at each probed distance
    recip_dist : array
        distance values between reciprocal pairs
    """
    unordered_probes = {}
    unordered_dist = {}
    # get all pairs from the two types which will have keys that are reciprocals of one another 
    # and values =  where there is a synapse
    pairs = probes_1 + probes_2
    probes = {(pair.pre_cell_id, pair.post_cell_id): pair.has_synapse for pair in pairs if pair.has_synapse is not None}
    probes_dist = {(pair.pre_cell_id, pair.post_cell_id): pair.distance for pair in pairs if pair.has_synapse is not None}

    # unorder the pairs above and count the number of synapses between them
    # 0 = no synapses
    # 1 = one-way synapse
    # 2 = reciprocal synapse
    for cell_ids,connected in probes.items():
        cell_ids = tuple(sorted(cell_ids))
        unordered_probes[cell_ids] = int(connected) + unordered_probes.get(cell_ids, 0)
        unordered_dist[cell_ids] = probes_dist.get(cell_ids, None)
        
    # set to True unordered pairs that have 2 synapses
    unordered_recip_connections = [x==2 for x in unordered_probes.values()]
    recip_conn = np.asarray(unordered_recip_connections, dtype='bool')
    recip_dist = np.asarray([x for x in unordered_dist.values()], dtype='float')
    
    connected1 = np.asarray([pair.has_synapse for pair in probes_1], dtype='bool')
    distance1 = np.asarray([pair.distance for pair in probes_1], dtype='float')
    _, cp_1, lower_1, upper_1 = connectivity_profile(connected1, distance1, bin_edges)

    connected2 = np.asarray([pair.has_synapse for pair in probes_2], dtype='bool')
    distance2 = np.asarray([pair.distance for pair in probes_2], dtype='float')
    _, cp_2, lower_2, upper_2 = connectivity_profile(connected2, distance2, bin_edges)

    _, cp_r, lower_r, upper_r = connectivity_profile(recip_conn, recip_dist, bin_edges)
    
    norm_cp_r = cp_r / (cp_1 * cp_2)

    return norm_cp_r, recip_conn, recip_dist


class CorrectionMetricFunctions:
    def pre_axon_length(pair):
        cell_morph = pair.pre_cell.morphology
        if cell_morph is None:
            return np.nan
        if cell_morph.axon_trunc_distance is not None:
            return cell_morph.axon_trunc_distance
        elif cell_morph.axon_truncation == 'intact':
            return 200e-6

    def avg_pair_depth(pair):
        if pair.pre_cell.depth is None or pair.post_cell.depth is None:
            return np.nan
        avg_depth = np.mean([pair.pre_cell.depth, pair.post_cell.depth])
        if avg_depth < 0 or avg_depth > 300e-6:
            return np.nan
        return avg_depth

    def n_test_spikes(pair):
        syn_type = pair.pre_cell.cell_class_nonsynaptic
        if syn_type not in ['ex', 'in']:
            return np.nan

        p = db.query(db.Pair).filter(db.Pair.id==pair.id).all()[0]
        n_spikes = getattr(p, 'n_%s_test_spikes' % syn_type)
        n_spikes = n_spikes if n_spikes <= 800 else 800
        return n_spikes            

    def baseline_rms_noise(pair):
        post_cell = pair.post_cell
        q = db.query(db.PatchClampRecording)
        q = q.join(db.Recording).join(db.Electrode).join(db.Cell)
        q = q.filter(db.Cell.id==post_cell.id).filter(db.PatchClampRecording.clamp_mode=='ic')
        pcrs = q.all()
        pcr_noise = [pcr.baseline_rms_noise for pcr in pcrs if pcr.qc_pass]
        if len(pcr_noise) > 1:
            return np.mean(pcr_noise)
        else:
            return np.nan

    def detection_power(pair):
        n_spikes = CorrectionMetricFunctions.n_test_spikes(pair)
        baseline_noise = CorrectionMetricFunctions.baseline_rms_noise(pair)
        if n_spikes != np.nan and baseline_noise != np.nan:
            dp = np.log10(np.sqrt(n_spikes) / baseline_noise)
            if np.isfinite(dp):
                return dp
            else:
                return np.nan
        else:
            return np.nan


def ei_correct_connectivity(ei_classes, correction_metrics, pairs):
    # correction_metrics is a dictionary where the key is the name of the metric which corresponds to either a column in the 
    # SynPhys DB or is a callable function and the values are a dictionary identifying the fit model and any associated parameters
    
    correction_fits = {}
    for class_name, cell_class in ei_classes.items():

        class_pairs = [p for p in pairs if p.pre_cell.cell_class_nonsynaptic==cell_class.name]

        probed_pairs = [p for p in class_pairs if pair_was_probed(p, cell_class.output_synapse_type)]
        connections = np.array([p.has_synapse for p in probed_pairs], dtype=bool)

        correction_fits[class_name] = []
        for metric, opts in correction_metrics.items():
            pair_metric = np.array([getattr(p, metric) for p in probed_pairs], dtype=float)    
            mask = np.isfinite(pair_metric) & np.isfinite(connections)
            model = opts['model']
            model_opts = {k:v for k, v in opts.items() if k != 'model'}

            metric_fit = model.fit(pair_metric[mask], connections[mask], **model_opts)
            correction_fits[class_name].append(metric_fit)
    
    correction_parameters = [[fit.fit_result.x for fit in fits] for fits in correction_fits.values()]
    correction_functions = [metric['model'].correction_func for metric in correction_metrics.values()]
    corr_model = CorrectionModel(0.1, correction_metrics.keys(), correction_functions, correction_parameters)
    
    return corr_model
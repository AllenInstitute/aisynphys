# coding: utf8
from __future__ import print_function, division
import numpy as np
import numba
import scipy.optimize


# lets us quickly disable jit for debugging:
def _fake_jit(**kwds):
    return lambda fn: fn
#jit = _fake_jit
jit = numba.jit


class StochasticReleaseModel(object):
    """A model of stochastic synaptic release used for determining optimal model parameters that describe
    empirically measured synaptic responses.
    
    Synaptic strength changes moment to moment based on the prior history of action potentials at 
    the presynaptic terminal. Many models have been published previously that attempt to capture this relationship.
    However, synaptic strength also depends on the synapse's prior history of vesicle release. This model uses
    both spike timing and response amplitude to compare a series of evoked synaptic events against the distribution
    of likely amplitudes predicted by the model.
    
    Usually, synaptic response data alone is not sufficient to fully constrain all parameters in a release model.
    This model is intended to be used to search large paremeter spaces to determine the subspace of parameters
    that are consistent with the measured data.
    
    Parameters
    ----------
    params : dict
        A dictionary of parameters specifying the behavior of the model:

        - n_release_sites (int) : Number of synaptic release zones
        - base_release_probability (float) : Resting-state synaptic release probability (0.0-1.0)
        - mini_amplitude (float) : Mean PSP amplitude evoked by a single vesicle release
        - mini_amplitude_cv (float) : Coefficient of variation of PSP amplitude evoked from single vesicle releases
        - depression_amount (float) : Amount of depression (0.0-1.0) to apply per spike. The special value -1 enables vesicle depletion rather than Pr depression.
        - depression_tau (float) : Time constant for recovery from depression or vesicle depletion
        - facilitation_amount (float) : Release probability facilitation per spike (0.0-1.0)
        - facilitation_tau (float) : Time constant for facilitated release probability to recover toward resting state
        - measurement_stdev (float) : Extra variance in PSP amplitudes purely as a result of membrane noise / measurement error
    """
    
    result_dtype = [
        ('spike_time', float),
        ('amplitude', float),
        ('expected_amplitude', float),
        ('likelihood', float),
    ]
    
    state_dtype = [
        ('vesicle_pool', float),
        ('depression', float),
        ('facilitation', float),
        ('release_probability', float),
    ]

    param_names = [
        'n_release_sites',
        'base_release_probability',
        'mini_amplitude',
        'mini_amplitude_cv',
        'depression_amount',
        'depression_tau',
        'facilitation_amount',
        'facilitation_tau',
        'measurement_stdev',
    ]
        
    def __init__(self, params):
        for k in params:
            if k not in self.param_names:
                raise ValueError("Unknown parameter name %r" % k)
        self.params = params

        # How long to wait after a NaN event before the model begins accumulating likelihood values again
        self.missing_event_penalty = 0.0

    def optimize_mini_amplitude(self, spike_times, amplitudes, event_meta=None, show=False):
        """Given a set of spike times and amplitudes, optimize the mini_amplitude parameter
        to produce the highest likelihood model.
        
        Returns the output of run_model for the best model.
        """
        params = self.params.copy()
        
        init_amp = estimate_mini_amplitude(amplitudes, params)
        params['mini_amplitude'] = init_amp
        init_result = self.run_model(spike_times, amplitudes, params, event_meta=event_meta)
        mean_amp = np.nanmean(amplitudes)
        if show:
            print("========== Optimize mini amplitude ==============")
            for k,v in params.items():
                print("   %s: %s" % (k, v))
            print("   initial guess:", init_amp)
            print("   mean amp:", mean_amp)
        ratio = mean_amp / np.nanmean(init_result.result['expected_amplitude'])
        init_amp *= ratio
        if show:
            print("   corrected amp 1:", init_amp)
        
        # in some cases, init_amp is way too large; force these back down:
        init_amp = min(init_amp, mean_amp) if mean_amp > 0 else max(init_amp, mean_amp)
        if show:
            print("   corrected amp 2:", init_amp)
             
        # self.params['mini_amplitude'] = params['mini_amplitude']
        result = self.optimize(spike_times, amplitudes, optimize={'mini_amplitude': (init_amp, init_amp*0.01, init_amp*100)}, event_meta=event_meta)
        if show:
            print("   optimized amp:", result.optimized_params['mini_amplitude'])
        result.optimization_info = {'init_amp': init_amp / ratio, 'ratio': ratio, 'corrected_amp': init_amp, 'init_likelihood': init_result.likelihood}
        return result
    
    def optimize(self, spike_times, amplitudes, optimize, event_meta=None):
        """Optimize specific parameters to maximize the model likelihood.

        This method updates the attributes for any optimized parameters and returns the
        best result from run_model().
        
        Parameters
        ----------
        spike_times : array
            Times (in seconds) of presynaptic spikes in ascending order
        amplitudes : array
            Evoked PSP/PSC amplitudes for each spike listed in *spike_times*. Amplitudes may be
            NaN to indicate that the event should be ignored (usually because the spike could not be
            detected or the amplitude could not be measured). Any events within 10 seconds following
            a NaN will update the model as usual, but will not be included in the likelihood estimate.
        optimize : dict | callable
            A dictionary of {'parameter_name': (init, min, max)} that specifies model
            parameters to be optimized. Alternatively, this may be a function that accepts 
            a dictionary of fixed model parameters and returns the dictionary of parameters
            to be optimized, as described above.
        """
        params = self.params.copy()
        if callable(optimize):
            optimize = optimize(self, params)

        init = [val[0] for val in optimize.values()]
        bounds = [sorted(val[1:3]) for val in optimize.values()]
        results = {}
        
        def fn(x):
            opt_params = params.copy()
            for i,k in enumerate(optimize.keys()):
                opt_params[k] = np.clip(x[i], *bounds[i])
            res = self.run_model(spike_times, amplitudes, opt_params, event_meta=event_meta)
            results[tuple(x)] = res
            return -res.likelihood
        
        best = scipy.optimize.minimize(fn, x0=init, 
            method="Nelder-Mead", options={'fatol': 0.01}  # no bounds, can't set initial step?
            #method='BFGS', bounds=bounds, options={'gtol': 1, 'eps': 10e-6}  # oscillates
            #method='Powell', options={'ftol': 0.01}  # not efficient; can't set initial step
            #method='CG', options={'eps': 10e-6}  # not efficient
            #method='Newton-CG',  # requires jacobian
            #method='L-BFGS-B'  # not efficient
            #method='COBYLA', options={'rhobeg': -100e-6}  # fast but oscillates, misses peak
        )
        
        best_result = results[tuple(best.x.flatten())]

        # take optimized params out of result.params and put them in result.optimized_params instead
        # (to ensure we can re-run in the same way)
        best_result.params = self.params.copy()
        best_result.optimized_params = {k:best.x[i] for i,k in enumerate(optimize.keys())}
        best_result.optimization_init = optimize
        best_result.optimization_result = best
        best_result.optimization_path = {
            'mini_amplitude': [k[0] for k in results.keys()], 
            'likelihood': [v.likelihood for v in results.values()]
        }
        
        # update attributes with best result
        self.params = self.params.copy()
        for i,k in enumerate(optimize.keys()):
            self.params[k] = best.x[i]
                    
        return best_result
    
    def run_model(self, spike_times, amplitudes, params=None, event_meta=None):
        """Run the stochastic release model with a specific set of spike times.
        
        This can be used two different ways: (1) compute a measure of the likelihood that *times* and *amplitudes*
        could be generated by a synapse with the current dynamic parameters, or (2) run the model to predict
        the behavior of a synapse given a set of presynaptic spike times.
        
        Parameters
        ----------
        spike_times : array
            Times (in seconds) of presynaptic spikes in ascending order
        amplitudes : array | str
            If this argument is an array, then it specifies the evoked PSP/PSC amplitudes for each spike 
            listed in *spike_times*. Amplitudes may be NaN to indicate that the event should be ignored 
            (usually because the spike could not be detected or the amplitude could not be measured). 
            Any events within 10 seconds following a NaN will update the model as usual, but will not 
            be included in the likelihood estimate.
            If this argument is a string, then the model will generate simulated amplitudes as it runs.
            The string may be "random" to randomly select amplitudes from the model distribution, or
            "expected" to use the expectation value of the predicted distribution at each spike time.
        params :  dict
            Dictionary of model parameter values to use during this run. By default, parameters
            are taken from self.params.
        event_meta : array
            Extra per-event metadata to be included in the model results
        
        Returns
        -------
        result : StochasticReleaseModelResult
        """
        if params is None:
            params = self.params
        use_expectation = False
        if isinstance(amplitudes, str):
            assert amplitudes in ('expected', 'random'), "amplitudes argument must be ndarray, 'expected', or 'random'"
            use_expectation = amplitudes == 'expected'
            amplitudes = np.array([])
        
        assert params['n_release_sites'] < 67, "For n_release_sites > 66 we need to use scipy.special.binom instead of the optimized binom_coeff"

        result = np.empty(len(spike_times), dtype=self.result_dtype)
        pre_spike_state = np.full(len(spike_times), np.nan, dtype=self.state_dtype)
        post_spike_state = np.full(len(spike_times), np.nan, dtype=self.state_dtype)
        
        try:
            self._run_model(
                spike_times=spike_times, 
                amplitudes=amplitudes, 
                result=result, 
                pre_spike_state=pre_spike_state, 
                post_spike_state=post_spike_state, 
                missing_event_penalty=self.missing_event_penalty,
                use_expectation=use_expectation,
                **params,
            )
        except Exception as exc:
            print("Error with params:", params)
            raise
        
        # scalar representation of overall likelihood
        likelihood = np.exp(np.nanmean(np.log(result['likelihood'] + 0.1)))
        
        return StochasticReleaseModelResult(
            result=result, 
            pre_spike_state=pre_spike_state,
            post_spike_state=post_spike_state,
            likelihood=likelihood,
            params=params,
            model=self,
            event_meta=event_meta,
        )
            
    def likelihood(self, amplitudes, state, params=None):
        """Estimate the probability density of seeing a particular *amplitude*
        given a number of *available_vesicles*.
        """
        if params is None:
            params = self.params.copy()
        
        available_vesicles = int(np.clip(np.round(state['vesicle_pool']), 0, params['n_release_sites']))
        return release_likelihood(amplitudes, available_vesicles, state['release_probability'], params['mini_amplitude'], params['mini_amplitude_cv'], params['measurement_stdev'])

    @staticmethod
    @jit(nopython=True)
    def _run_model( spike_times, 
                    amplitudes,
                    result,
                    pre_spike_state,
                    post_spike_state, 
                    missing_event_penalty,
                    use_expectation,
                    n_release_sites,
                    base_release_probability,
                    mini_amplitude,
                    mini_amplitude_cv,
                    depression_amount,
                    depression_tau,
                    facilitation_amount,
                    facilitation_tau,
                    measurement_stdev,
                    ):

        have_amps = len(amplitudes) > 0

        # initialize state parameters:
        # vesicle_pool is a float as a means of avoiding the need to model stochastic vesicle docking;
        # we just assume that recovery is a continuous process. (maybe we should just call this "neurotransmitter"
        # instead)
        vesicle_pool = n_release_sites
        facilitation = 0.0
        depression = 0.0
        
        # vesicle depletion is enabled if depression_amount is -1
        use_vesicle_depletion = depression_amount == -1

        previous_t = spike_times[0]
        last_nan_time = -np.inf

        for i,t in enumerate(spike_times):
            if have_amps:
                amplitude = amplitudes[i]
            
            # nan amplitude means that a stimulus occurred here but it is uncertain whether a spike
            # was evoked. In this case, we skip over the event and potentially incur a timeout before
            # using the likelihood from future events.
            if have_amps and np.isnan(amplitude):
                expected_amplitude = np.nan
                likelihood = np.nan
                last_nan_time = t

            else:
                dt = t - previous_t
                previous_t = t

                # recover vesicles up to the current timepoint if vesicle depletion is enabled
                if use_vesicle_depletion:
                    v_recovery = np.exp(-dt / depression_tau)
                    vesicle_pool += (n_release_sites - vesicle_pool) * (1.0 - v_recovery)

                    # bounded, integer number of available vesicles used for binomial
                    effective_available_vesicles = max(0, min(n_release_sites, int(np.round(vesicle_pool))))
                    # unbounded, continuous number of available vesicles used for calculating expectation value
                    effective_available_vesicle = max(0, vesicle_pool)
                else:
                    effective_available_vesicles = int(vesicle_pool)
                    effective_available_vesicle = vesicle_pool

                    # recover depression mechanism only if vesicle depletion is disabled
                    depression *= np.exp(-dt / depression_tau)
                
                # recover facilitation mechanism
                facilitation *= np.exp(-dt / facilitation_tau)

                # calculate release probability
                release_probability = (1 - depression) * (base_release_probability + (1 - base_release_probability) * facilitation)

                # predict most likely amplitude for this spike (just for show)
                expected_amplitude = release_expectation_value(
                    effective_available_vesicle,
                    release_probability,
                    mini_amplitude,
                )
                
                # Generate response amplitudes if they were not already supplied
                if not have_amps:
                    if use_expectation:
                        # run model forward based on expectation value
                        amplitude = expected_amplitude
                    else:
                        # select a random amplitude from distribution
                        amplitude = release_random_value(
                            effective_available_vesicles,
                            release_probability, 
                            mini_amplitude,
                            mini_amplitude_cv,
                            measurement_stdev
                        )


                # Measure likelihood
                if t - last_nan_time < missing_event_penalty:
                    # ignore likelihood for this event if it was too close to an unmeasurable response
                    likelihood = np.nan
                else:
                    # measure likelihood of seeing this response amplitude
                    likelihood = release_likelihood_scalar(amplitude, effective_available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev)
                # prof('likelihood')


                # record model state immediately before spike
                pre_spike_state[i]['vesicle_pool'] = vesicle_pool
                pre_spike_state[i]['release_probability'] = release_probability
                pre_spike_state[i]['depression'] = depression
                pre_spike_state[i]['facilitation'] = facilitation

                # release vesicles
                if use_vesicle_depletion:
                    # note: we allow vesicle_pool to become negative because this helps to ensure
                    # that the overall likelihood will be low for such models; however, the value
                    # used for binomial is bounded at 0
                    depleted_vesicle = amplitude / mini_amplitude
                    vesicle_pool -= depleted_vesicle
                    assert np.isfinite(vesicle_pool)
                else:
                    # apply spike-induced depression only if vesicle depletion is disabled
                    depression += (1 - depression) * depression_amount

                # apply spike-induced facilitation
                facilitation += (1 - facilitation) * facilitation_amount

                # prof('update state')

                # record model state immediately after spike
                post_spike_state[i]['vesicle_pool'] = vesicle_pool
                post_spike_state[i]['release_probability'] = release_probability
                post_spike_state[i]['depression'] = depression
                post_spike_state[i]['facilitation'] = facilitation
                # prof('record')
            
            # record results
            result[i]['spike_time'] = t
            result[i]['amplitude'] = amplitude
            result[i]['expected_amplitude'] = expected_amplitude
            result[i]['likelihood'] = likelihood


class StochasticReleaseModelResult:
    """Contains the results of StochasticReleaseModel.run_model in several attributes:

    Attributes
    ----------
    result : array
        Array of model results, one record per input spike.
        Contains fields: spike_time, amplitude, expected_amplitude, likelihood
    pre_spike_state : array
        State variable values immediately before each spike
    post_spike_state : array
        State variable values immediately after each spike
    likelihood : float
        The estimated likelihood for this model run
    params : dict
        A dictionary of parameters used to generate the model result
    model : StochasticReleaseModel
        The model instance that generated this result
    event_meta : array
        Per-event metadata, mainly regarding stimulus structure
    optimized_params : dict | None
        Parameters that were optimized by running the model
    """
    pickle_attributes = ['likelihood', 'params', 'optimized_params']

    def __init__(self, result, pre_spike_state, post_spike_state, likelihood, params, model, optimized_params=None, event_meta=None):
        self.result = result
        self.pre_spike_state = pre_spike_state
        self.post_spike_state = post_spike_state
        self.likelihood = likelihood
        self.params = params
        self.model = model
        self.event_meta = event_meta

        # filled in later by optimization routine
        self.optimized_params = optimized_params or {}

    @property
    def all_params(self):
        """A dictionary containing the combination of self.params and self.optimized_params; can be used 
        to re-run the model with the result of the parameter optimization fixed.
        """
        p = self.params.copy()
        p.update(self.optimized_params)
        return p

    def model_amp_distribution(self, amp_values, index):
        """Return the model exepcted amplitude probability distribution immediately prior to *index*.

        If the model state is not well defined at *index*, then return None.
        """
        state = self.pre_spike_state[index]
        if not np.all(np.isfinite(tuple(state))):
            return None
        return self.model.likelihood(amp_values, state)

    def avg_model_amp_distribution(self, amp_values=None, index=None):
        """Return the average model exepcted amplitude probability distribution.

        Parameters
        ----------
        index : array | None
            Event indices at which the model distribution should be sampled for averaging.
            If None, then all events are included in the average.
        amp_values : array
            Amplitude values at which the model distribution will be evaluated.
        """
        if index is None:
            index = range(len(self.pre_spike_state))

        dist = None
        count = 0
        for i in index:
            d = self.model_amp_distribution(amp_values, i)
            if d is None:
                continue
            if dist is None:
                dist = d
            else:
                dist += d
                count += 1

        return dist / count

    def events_by_recording(self, require_contiguous=True):
        """Return event ids grouped by recording.

        Parameters
        ----------
        require_contiguous : bool
            If True, only return the IDs of events for which a complete train of events is available up to that point.
            (for example, if we have a 12-pulse train but pule 5 is missing, then only pulses 1-4 will be included in the output)

        Returns
        -------
        recs : list
            Each item in the returned list is a list of event IDs, all of which originated from the same recording.
        """
        spikes = self.result['spike_time']
        amps = self.result['amplitude']
        recs = []
        last_rec_id = None
        last_pulse_n = 0
        for i in range(len(self.result)):
            rec_id = self.event_meta['sync_rec_ext_id'][i]
            if rec_id != last_rec_id:
                last_rec_id = rec_id
                if len(recs) == 0 or len(recs[-1]) > 0:
                    recs.append([])
                last_pulse_n = 0
            pulse_n = self.event_meta['pulse_number'][i]
            if not require_contiguous or (pulse_n == last_pulse_n + 1 and np.isfinite(spikes[i]) and np.isfinite(amps[i])):
                recs[-1].append(i)
                last_pulse_n = pulse_n
        if len(recs[-1]) == 0:
            recs.pop(-1)
        return recs

    def events_by_stimulus(self, require_contiguous=True):
        """Return a structure giving event IDs grouped by stimulus parameters.

        Structure like::

            {induction_freq: {recovery_delay: [[a1, a2, ..a12], [b1, b2, ..b12], ...], ...}, ...}
        """
        recs = self.events_by_recording(require_contiguous)
        meta = self.event_meta
        
        trains = {}
        for rec in recs:
            i = rec[0]
            ind_f = meta['induction_frequency'][i]
            rec_d = meta['recovery_delay'][i]
            ind_trains = trains.setdefault(ind_f, {})
            rec_trains = ind_trains.setdefault(rec_d, [])
            rec_trains.append(rec)

        return trains

    def __getstate__(self):
        return {k: getattr(self, k) for k in self.pickle_attributes}

    def __setstate__(self, state):
        for k,v in state.items():
            setattr(self, k, v)


@jit(nopython=True)
def release_likelihood(amplitudes, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev):
    """Return a measure of the likelihood that a synaptic response will have certain amplitude(s),
    given the state parameters for the synapse.
    
    Parameters
    ----------
    amplitudes : array
        The amplitudes for which likelihood values will be returned
    available_vesicles : int
        Number of vesicles available for release
    release_probability : float
        Probability for each available vesicle to be released
    mini_amplitude : float
        Mean amplitude of response evoked by a single vesicle release
    mini_amplitude_cv : float
        Coefficient of variation of response amplitudes evoked by a single vesicle release
    measurement_stdev : float
        Standard deviation of response amplitude measurement errors
        
        
    For each value in *amplitudes*, we calculate the likelihood that a synapse would evoke a response
    of that amplitude. Likelihood is calculated as the probability denstity of a gaussian mixture, as follows:
    
    1. Given the number of vesicles available to be released (nV) and the release probability (pR), determine
       the probability that each possible number of vesicles (nR) will be released using the binomial distribution
       probability mass function. For example, if there are 3 vesicles available and the release probability is
       0.1, then the possibilities are:
           vesicles released (nR)    probability
                                0    0.729
                                1    0.243
                                2    0.27
                                3    0.001
    2. For each possible number of released vesicles, calculate the likelihood that this possibility could
       evoke a response of the tested amplitude. This is calculated using the Gaussian probability distribution 
       function where µ = nR * mini_amplitude and σ = sqrt((mini_amplitude * mini_amplitude_cv)^2 * nR + measurement_stdev)
    3. The total likelihood is the sum of likelihoods for all possible values of nR.
    """
    return np.array([
        release_likelihood_scalar(amplitude, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev)
        for amplitude in amplitudes])


@jit(nopython=True)
def release_likelihood_scalar(amplitude, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev):
    """Same as release_likelihood, but optimized for a scalar amplitude argument"""
    n_vesicles = np.arange(available_vesicles + 1)
    
    # probability of releasing n_vesicles given available_vesicles and release_probability
    p_n = binom_pmf(available_vesicles, release_probability, n_vesicles)
    
    # expected amplitude for n_vesicles
    amp_mean = n_vesicles * mini_amplitude
    
    # amplitude stdev increases by sqrt(n) with number of released vesicles
    amp_stdev = ((mini_amplitude * mini_amplitude_cv)**2 * n_vesicles + measurement_stdev**2) ** 0.5
    
    # distributions of amplitudes expected for n_vesicles
    amp_prob = p_n * normal_pdf(amp_mean, amp_stdev, amplitude)
    
    # sum all distributions across n_vesicles
    likelihood = amp_prob.sum()
    
    assert likelihood >= 0
    return likelihood


# def release_distribution(available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev):
#     """Return the release amplitude distribution defined by the arguments.
#     """
#     # calculate positive and negate at the end if needed.
#     sign = mini_amplitude / abs(mini_amplitude)
#     mini_amplitude = abs(mini_amplitude)
    
#     mn = -measurement_stdev * 3
#     mx = max(0, available_vesicles) * mini_amplitude + measurement_stdev * 3
#     n_samp = int(max(0, available_vesicles) + 1) * 20
#     amplitudes = np.linspace(mn, mx, n_samp)
#     da = amplitudes[1] - amplitudes[0]
#     return amplitudes * sign, release_likelihood(amplitudes, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev) * da


@jit(nopython=True)
def release_expectation_value(available_vesicles, release_probability, mini_amplitude):
    """Return the expectation value for the release amplitude distribution defined by the arguments.
    """
    return binom_mean(available_vesicles, release_probability) * mini_amplitude


@jit(nopython=True)
def release_random_value(
    available_vesicles,
    release_probability, 
    mini_amplitude,
    mini_amplitude_cv,
    measurement_stdev):
        n_vesicles = np.random.binomial(n=available_vesicles, p=release_probability)
        amp_mean = n_vesicles * mini_amplitude
        amp_stdev = ((mini_amplitude * mini_amplitude_cv)**2 * n_vesicles + measurement_stdev**2) ** 0.5
        return np.random.normal(loc=amp_mean, scale=amp_stdev)


@jit(nopython=True)
def normal_pdf(mu, sigma, x):
    """Probability density function of normal distribution
    
    Same as scipy.stats.norm(mu, sigma).pdf(x)
    """
    return (1.0 / (2 * np.pi * sigma**2))**0.5 * np.exp(- (x-mu)**2 / (2 * sigma**2))

#@functools.lru_cache(maxsize=2**14)
@jit(nopython=True)
def binom_pmf(n, p, k):
    """Probability mass function of binomial distribution
    
    Same as scipy.stats.binom(n, p).pmf(k), but much faster. Given *n* independent experiments,
    each with probability *p* of succeeding, return the probability of seeing *k* successful
    events (for each value given in the array *k*).

    Parameters
    ----------
    n : int
        Number of independent experiments
    p : float
        Probability of success per experiment
    k : array of int
        Numbers of successful exeriments for which to return the probability mass function
    """
    bc = np.array([binom_coeff(n, k1) for k1 in k])
    return bc * (p**k * (1 - p)**(n - k))


_binom_coeff_cache = np.fromfunction(scipy.special.binom, (67, 67)).astype(int)

@jit(nopython=True)
def binom_coeff(n, k):
    """Binomial coefficient: n! / (k! (n-k)!)
    
    Same as scipy.special.binom, but much faster and limited to n < 67.
    """
    # note: one cold imagine writing an optimized binomial coefficient function that
    # is not limited to n < 67, but having to look out for integer overflows slows us down.
    return _binom_coeff_cache[n, k]


@jit(nopython=True)
def binom_mean(n, p):
    """Expectation value of binomial distribution
    
    Same as stats.binom(n, p).mean(), but much-much faster.
    """
    return n * p


def estimate_mini_amplitude(amplitudes, params):
    avg_amplitude = np.nanmean(amplitudes)
    expected = release_expectation_value(params['n_release_sites'], params['base_release_probability'], mini_amplitude=1)
    init_amp = avg_amplitude / expected
    
    # takes care of cases where the first peak in the distribution is larger than any measured events;
    # otherwise this case is very difficult to optimize
    while abs(init_amp) > abs(avg_amplitude):
        init_amp /= 2

    return init_amp




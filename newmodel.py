import brian2 as b2
from brian2 import NeuronGroup, Synapses, PoissonInput, PoissonGroup, network_operation
from brian2.monitors import StateMonitor, SpikeMonitor, PopulationRateMonitor
from random import sample
import numpy.random as rnd
from neurodynex.tools import plot_tools
import numpy
import matplotlib.pyplot as plt
from math import floor
import time

b2.defaultclock.dt = 0.10 * b2.ms
def sim_decision_making_network_full(N_Excit = 384, N_Inhib = 96, weight_scaling_factor = 5.33,
                                t_stimulus_start = 100 * b2.ms, t_stimulus_duration = 9999 * b2.ms, coherence_level = 50.0,
                                stimulus_update_interval =30 * b2.ms, mu0_mean_stimulus_Hz = 25.,
                                stimulus_std_Hz = 7.,
                                N_extern = 1000, firing_rate_extern = 9.8 * b2.Hz,
                                w_pos = 1.70, f_Subpop_size = 0.25,  # .15 in publication [1]
                                max_sim_time = 2000. * b2.ms, stop_condition_rate = 17 * b2.Hz,
                                monitored_subset_size = 512):
    """

    Args:
        N_Excit (int): total number of neurons in the excitatory population
        N_Inhib (int): nr of neurons in the inhibitory populations
        weight_scaling_factor: When increasing the number of neurons by 2, the weights should be scaled down by 1/2
        t_stimulus_start (Quantity): time when the stimulation starts
        t_stimulus_duration (Quantity): duration of the stimulation
        coherence_level (int): coherence of the stimulus.
            Difference in mean between the PoissonGroups "left" stimulus and "right" stimulus
        stimulus_update_interval (Quantity): the mean of the stimulating PoissonGroups is
            re-sampled at this interval
        mu0_mean_stimulus_Hz (float): maximum mean firing rate of the stimulus if c=+1 or c=-1. Each neuron
            in the populations "Left" and "Right" receives an independent poisson input.
        stimulus_std_Hz (float): std deviation of the stimulating PoissonGroups.
        N_extern (int): nr of neurons in the stimulus independent poisson background population
        firing_rate_extern (int): firing rate of the stimulus independent poisson background population
        w_pos (float): Scaling (strengthening) of the recurrent weights within the
            subpopulations "Left" and "Right"
        f_Subpop_size (float): fraction of the neurons in the subpopulations "Left" and "Right".
            #left = #right = int(f_Subpop_size*N_Excit).
        max_sim_time (Quantity): simulated time.
        stop_condition_rate (Quantity): An optional stopping criteria: If not None, the simulation stops if the
            firing rate of either subpopulation "Left" or "Right" is above stop_condition_rate.
        monitored_subset_size (int): max nr of neurons for which a state monitor is registered.

    Returns:

        A dictionary with the following keys (strings):
        "rate_monitor_A", "spike_monitor_A", "voltage_monitor_A", "idx_monitored_neurons_A", "rate_monitor_B",
         "spike_monitor_B", "voltage_monitor_B", "idx_monitored_neurons_B", "rate_monitor_Z", "spike_monitor_Z",
         "voltage_monitor_Z", "idx_monitored_neurons_Z", "rate_monitor_inhib", "spike_monitor_inhib",
         "voltage_monitor_inhib", "idx_monitored_neurons_inhib"

    """
    startting_time = time.time()
    t_stimulus_end = t_stimulus_start + t_stimulus_duration

    N_Group_A = int(N_Excit * f_Subpop_size)  # size of the excitatory subpopulation sensitive to stimulus A
    N_Group_B = N_Group_A  # size of the excitatory subpopulation sensitive to stimulus B
    N_Group_Z = N_Excit - N_Group_A - N_Group_B  # (1-2f)Ne excitatory neurons do not respond to either stimulus.

    Cm_excit = 0.5 * b2.nF  # membrane capacitance of excitatory neurons
    G_leak_excit = 25.0 * b2.nS  # leak conductance
    E_leak_excit = -70.0 * b2.mV  # reversal potential
    v_spike_thr_excit = -50.0 * b2.mV  # spike condition
    v_reset_excit = -55.0 * b2.mV  # reset voltage after spike
    t_abs_refract_excit = 2.0 * b2.ms  # absolute refractory period

    # specify the inhibitory interneurons:
    # N_Inhib = 200
    Cm_inhib = 0.2 * b2.nF
    G_leak_inhib = 20.0 * b2.nS
    E_leak_inhib = -70.0 * b2.mV
    v_spike_thr_inhib = -50.0 * b2.mV
    v_reset_inhib = -55.0 * b2.mV
    t_abs_refract_inhib = 1.0 * b2.ms

    # specify the AMPA synapses
    E_AMPA = 0.0 * b2.mV
    tau_AMPA = 2.0 * b2.ms

    # specify the GABA synapses
    E_GABA = -70.0 * b2.mV
    tau_GABA = 5.0 * b2.ms

    # specify the NMDA synapses
    E_NMDA = 0.0 * b2.mV
    tau_NMDA_s = 100.0 * b2.ms
    tau_NMDA_x = 2.0 * b2.ms
    alpha_NMDA = 0.5 * b2.kHz

    # projections from the external population
    g_AMPA_extern2inhib = 1.62 * b2.nS
    g_AMPA_extern2excit = 2.1 * b2.nS

    # projectsions from the inhibitory populations
    g_GABA_inhib2inhib = weight_scaling_factor * 1.25 * b2.nS
    g_GABA_inhib2excit = weight_scaling_factor * 1.60 * b2.nS

    # projections from the excitatory population
    g_AMPA_excit2excit = weight_scaling_factor * 0.012 * b2.nS
    g_AMPA_excit2inhib = weight_scaling_factor * 0.015 * b2.nS
    g_NMDA_excit2excit = weight_scaling_factor * 0.040 * b2.nS
    g_NMDA_excit2inhib = weight_scaling_factor * 0.045 * b2.nS  # stronger projection to inhib.

    # weights and "adjusted" weights.
    w_neg = 1. - f_Subpop_size * (w_pos - 1.) / (1. - f_Subpop_size)
    # We use the same postsyn AMPA and NMDA conductances. Adjust the weights coming from different sources:
    w_ext2inhib = g_AMPA_extern2inhib / g_AMPA_excit2inhib
    w_ext2excit = g_AMPA_extern2excit / g_AMPA_excit2excit
    # other weights are 1
    # print("w_neg={}, w_ext2inhib={}, w_ext2excit={}".format(w_neg, w_ext2inhib, w_ext2excit))

    # Define the inhibitory population
    # dynamics:
    inhib_lif_dynamics = """
        s_NMDA_total : 1  # the post synaptic sum of s. compare with s_NMDA_presyn
        dv/dt = (
        - G_leak_inhib * (v-E_leak_inhib)
        - g_AMPA_excit2inhib * s_AMPA * (v-E_AMPA)
        - g_GABA_inhib2inhib * s_GABA * (v-E_GABA)
        - g_NMDA_excit2inhib * s_NMDA_total * (v-E_NMDA)/(1.0+1.0*exp(-0.062*v/volt)/3.57)
        )/Cm_inhib : volt (unless refractory)
        ds_AMPA/dt = -s_AMPA/tau_AMPA : 1
        ds_GABA/dt = -s_GABA/tau_GABA : 1
    """

    inhib_pop = NeuronGroup(
        N_Inhib, model = inhib_lif_dynamics,
        threshold = "v>v_spike_thr_inhib", reset = "v=v_reset_inhib", refractory = t_abs_refract_inhib,
        method = "rk2")
    # initialize with random voltages:
    inhib_pop.v = rnd.uniform(v_spike_thr_inhib / b2.mV - 4., high = v_spike_thr_inhib / b2.mV - 1., size = N_Inhib) * b2.mV

    # Specify the excitatory population:
    # dynamics:
    excit_lif_dynamics = """
        s_NMDA_total : 1  # the post synaptic sum of s. compare with s_NMDA_presyn
        dv/dt = (
        - G_leak_excit * (v-E_leak_excit)
        - g_AMPA_excit2excit * s_AMPA * (v-E_AMPA)
        - g_GABA_inhib2excit * s_GABA * (v-E_GABA)
        - g_NMDA_excit2excit * s_NMDA_total * (v-E_NMDA)/(1.0+1.0*exp(-0.062*v/volt)/3.57)
        )/Cm_excit : volt (unless refractory)
        ds_AMPA/dt = -s_AMPA/tau_AMPA : 1
        ds_GABA/dt = -s_GABA/tau_GABA : 1
        ds_NMDA/dt = -s_NMDA/tau_NMDA_s + alpha_NMDA * x * (1-s_NMDA) : 1
        dx/dt = -x/tau_NMDA_x : 1
    """

    # define the three excitatory subpopulations.
    # A: subpop receiving stimulus A
    excit_pop_A = NeuronGroup(N_Group_A, model = excit_lif_dynamics,
                              threshold = "v > v_spike_thr_excit", reset="v = v_reset_excit",
                              refractory = t_abs_refract_excit, method = "rk2")
    excit_pop_A.v = rnd.uniform(E_leak_excit / b2.mV, high = E_leak_excit / b2.mV + 5., size = excit_pop_A.N) * b2.mV

    # B: subpop receiving stimulus B
    excit_pop_B = NeuronGroup(N_Group_B, model = excit_lif_dynamics, threshold = "v > v_spike_thr_excit",
                              reset="v = v_reset_excit", refractory=t_abs_refract_excit, method="rk2")
    excit_pop_B.v = rnd.uniform(E_leak_excit / b2.mV, high = E_leak_excit / b2.mV + 5., size = excit_pop_B.N) * b2.mV
    # Z: non-sensitive
    excit_pop_Z = NeuronGroup(N_Group_Z, model = excit_lif_dynamics,
                              threshold = "v>v_spike_thr_excit", reset = "v=v_reset_excit",
                              refractory = t_abs_refract_excit, method="rk2")
    excit_pop_Z.v = rnd.uniform(v_reset_excit / b2.mV, high = v_spike_thr_excit / b2.mV - 1., size = excit_pop_Z.N) * b2.mV

    # now define the connections:
    # projections FROM EXTERNAL POISSON GROUP: ####################################################
    poisson2Inhib = PoissonInput(target=inhib_pop, target_var="s_AMPA",
                                 N=N_extern, rate=firing_rate_extern, weight=w_ext2inhib)
    poisson2A = PoissonInput(target=excit_pop_A, target_var="s_AMPA",
                             N=N_extern, rate=firing_rate_extern, weight=w_ext2excit)

    poisson2B = PoissonInput(target=excit_pop_B, target_var="s_AMPA",
                             N=N_extern, rate=firing_rate_extern, weight=w_ext2excit)
    poisson2Z = PoissonInput(target=excit_pop_Z, target_var="s_AMPA",
                             N=N_extern, rate=firing_rate_extern, weight=w_ext2excit)

    ###############################################################################################

    # GABA projections FROM INHIBITORY population: ################################################
    syn_inhib2inhib = Synapses(inhib_pop, target=inhib_pop, on_pre="s_GABA += 1.0", delay=0.5 * b2.ms)
    syn_inhib2inhib.connect(p=1.)
    syn_inhib2A = Synapses(inhib_pop, target=excit_pop_A, on_pre="s_GABA += 1.0", delay=0.5 * b2.ms)
    syn_inhib2A.connect(p=1.)
    syn_inhib2B = Synapses(inhib_pop, target=excit_pop_B, on_pre="s_GABA += 1.0", delay=0.5 * b2.ms)
    syn_inhib2B.connect(p=1.)
    syn_inhib2Z = Synapses(inhib_pop, target=excit_pop_Z, on_pre="s_GABA += 1.0", delay=0.5 * b2.ms)
    syn_inhib2Z.connect(p=1.)
    ###############################################################################################

    # AMPA projections FROM EXCITATORY A: #########################################################
    syn_AMPA_A2A = Synapses(excit_pop_A, target=excit_pop_A, on_pre="s_AMPA += w_pos", delay=0.5 * b2.ms)
    syn_AMPA_A2A.connect(p=1.)
    syn_AMPA_A2B = Synapses(excit_pop_A, target=excit_pop_B, on_pre="s_AMPA += w_neg", delay=0.5 * b2.ms)
    syn_AMPA_A2B.connect(p=1.)
    syn_AMPA_A2Z = Synapses(excit_pop_A, target=excit_pop_Z, on_pre="s_AMPA += w_neg", delay=0.5 * b2.ms)
    syn_AMPA_A2Z.connect(p=1.)
    syn_AMPA_A2inhib = Synapses(excit_pop_A, target=inhib_pop, on_pre="s_AMPA += 1.0", delay=0.5 * b2.ms)
    syn_AMPA_A2inhib.connect(p=1.)
    ###############################################################################################

    # AMPA projections FROM EXCITATORY B: #########################################################
    syn_AMPA_B2A = Synapses(excit_pop_B, target=excit_pop_A, on_pre="s_AMPA += w_neg", delay=0.5 * b2.ms)
    syn_AMPA_B2A.connect(p=1.)
    syn_AMPA_B2B = Synapses(excit_pop_B, target=excit_pop_B, on_pre="s_AMPA += w_pos", delay=0.5 * b2.ms)
    syn_AMPA_B2B.connect(p=1.)
    syn_AMPA_B2Z = Synapses(excit_pop_B, target=excit_pop_Z, on_pre="s_AMPA += w_neg", delay=0.5 * b2.ms)
    syn_AMPA_B2Z.connect(p=1.)
    syn_AMPA_B2inhib = Synapses(excit_pop_B, target=inhib_pop, on_pre="s_AMPA += 1.0", delay=0.5 * b2.ms)
    syn_AMPA_B2inhib.connect(p=1.)
    ###############################################################################################

    # AMPA projections FROM EXCITATORY Z: #########################################################
    syn_AMPA_Z2A = Synapses(excit_pop_Z, target=excit_pop_A, on_pre="s_AMPA += w_neg", delay=0.5 * b2.ms)
    syn_AMPA_Z2A.connect(p=1.)
    syn_AMPA_Z2B = Synapses(excit_pop_Z, target=excit_pop_B, on_pre="s_AMPA += w_neg", delay=0.5 * b2.ms)
    syn_AMPA_Z2B.connect(p=1.)
    syn_AMPA_Z2Z = Synapses(excit_pop_Z, target=excit_pop_Z, on_pre="s_AMPA += 1.0", delay=0.5 * b2.ms)
    syn_AMPA_Z2Z.connect(p=1.)
    syn_AMPA_Z2inhib = Synapses(excit_pop_Z, target=inhib_pop, on_pre="s_AMPA += 1.0", delay=0.5 * b2.ms)
    syn_AMPA_Z2inhib.connect(p=1.)
    ###############################################################################################

    # NMDA projections FROM EXCITATORY to INHIB, A,B,Z
    @network_operation()
    def update_nmda_sum():
        sum_sNMDA_A = sum(excit_pop_A.s_NMDA)
        sum_sNMDA_B = sum(excit_pop_B.s_NMDA)
        sum_sNMDA_Z = sum(excit_pop_Z.s_NMDA)
        # note the _ at the end of s_NMDA_total_ disables unit checking
        inhib_pop.s_NMDA_total_ = (1.0 * sum_sNMDA_A + 1.0 * sum_sNMDA_B + 1.0 * sum_sNMDA_Z)
        excit_pop_A.s_NMDA_total_ = (w_pos * sum_sNMDA_A + w_neg * sum_sNMDA_B + w_neg * sum_sNMDA_Z)
        excit_pop_B.s_NMDA_total_ = (w_neg * sum_sNMDA_A + w_pos * sum_sNMDA_B + w_neg * sum_sNMDA_Z)
        excit_pop_Z.s_NMDA_total_ = (w_neg * sum_sNMDA_A + w_neg * sum_sNMDA_B + 1.0 * sum_sNMDA_Z)

    # set a self-recurrent synapse to introduce a delay when updating the intermediate
    # gating variable x
    syn_x_A2A = Synapses(excit_pop_A, excit_pop_A, on_pre="x += 1.", delay=0.5 * b2.ms)
    syn_x_A2A.connect(j="i")
    syn_x_B2B = Synapses(excit_pop_B, excit_pop_B, on_pre="x += 1.", delay=0.5 * b2.ms)
    syn_x_B2B.connect(j="i")
    syn_x_Z2Z = Synapses(excit_pop_Z, excit_pop_Z, on_pre="x += 1.", delay=0.5 * b2.ms)
    syn_x_Z2Z.connect(j="i")
    ###############################################################################################

    # Define the stimulus: two PoissonInput with time time-dependent mean.
    poissonStimulus2A = PoissonGroup(N_Group_A, 0. * b2.Hz)
    syn_Stim2A = Synapses(poissonStimulus2A, excit_pop_A, on_pre="s_AMPA+=w_ext2excit")
    syn_Stim2A.connect(j="i")
    poissonStimulus2B = PoissonGroup(N_Group_B, 0. * b2.Hz)
    syn_Stim2B = Synapses(poissonStimulus2B, excit_pop_B, on_pre="s_AMPA+=w_ext2excit")
    syn_Stim2B.connect(j="i")

    @network_operation(dt=stimulus_update_interval)
    def update_poisson_stimulus(t):
        if t >= t_stimulus_start and t < t_stimulus_end:
            offset_A = mu0_mean_stimulus_Hz * (0.5 + 0.5 * 0.01 * coherence_level)
            offset_B = mu0_mean_stimulus_Hz * (0.5 - 0.5 * 0.01 * coherence_level)

            rate_A = numpy.random.normal(offset_A, stimulus_std_Hz)
            rate_A = (max(0, rate_A)) * b2.Hz  # avoid negative rate
            rate_B = numpy.random.normal(offset_B, stimulus_std_Hz)
            rate_B = (max(0, rate_B)) * b2.Hz

            poissonStimulus2A.rates = rate_A
            poissonStimulus2B.rates = rate_B
            # print("stim on. rate_A= {}, rate_B = {}".format(rate_A, rate_B))
        else:
            # print("stim off")
            poissonStimulus2A.rates = 0.
            poissonStimulus2B.rates = 0.

    ###############################################################################################

    def get_monitors(pop, monitored_subset_size):
        monitored_subset_size = min(monitored_subset_size, pop.N)
        idx_monitored_neurons = sample(range(pop.N), monitored_subset_size)
        rate_monitor = PopulationRateMonitor(pop)
        # record parameter: record=idx_monitored_neurons is not supported???
        spike_monitor = SpikeMonitor(pop, record=idx_monitored_neurons)
        voltage_monitor = StateMonitor(pop, "v", record=idx_monitored_neurons)
        return rate_monitor, spike_monitor, voltage_monitor, idx_monitored_neurons

    # collect data of a subset of neurons:
    rate_monitor_inhib, spike_monitor_inhib, voltage_monitor_inhib, idx_monitored_neurons_inhib = \
        get_monitors(inhib_pop, monitored_subset_size)

    rate_monitor_A, spike_monitor_A, voltage_monitor_A, idx_monitored_neurons_A = \
        get_monitors(excit_pop_A, monitored_subset_size)

    rate_monitor_B, spike_monitor_B, voltage_monitor_B, idx_monitored_neurons_B = \
        get_monitors(excit_pop_B, monitored_subset_size)

    rate_monitor_Z, spike_monitor_Z, voltage_monitor_Z, idx_monitored_neurons_Z = \
        get_monitors(excit_pop_Z, monitored_subset_size)

    if stop_condition_rate is None:
        b2.run(max_sim_time)
    else:
        sim_sum = 0. * b2.ms
        sim_batch = 100. * b2.ms
        samples_in_batch = int(floor(sim_batch / b2.defaultclock.dt))
        avg_rate_in_batch = 0
        while (sim_sum < max_sim_time) and (avg_rate_in_batch < stop_condition_rate):
            b2.run(sim_batch)
            avg_A = numpy.mean(rate_monitor_A.rate[-samples_in_batch:])
            avg_B = numpy.mean(rate_monitor_B.rate[-samples_in_batch:])
            avg_rate_in_batch = max(avg_A, avg_B)
            sim_sum += sim_batch

    print("Time elapsed =", time.time() - startting_time, 'seconds')
    ret_vals = dict()

    ret_vals["rate_monitor_A"] = rate_monitor_A
    ret_vals["spike_monitor_A"] = spike_monitor_A
    ret_vals["voltage_monitor_A"] = voltage_monitor_A
    ret_vals["idx_monitored_neurons_A"] = idx_monitored_neurons_A

    ret_vals["rate_monitor_B"] = rate_monitor_B
    ret_vals["spike_monitor_B"] = spike_monitor_B
    ret_vals["voltage_monitor_B"] = voltage_monitor_B
    ret_vals["idx_monitored_neurons_B"] = idx_monitored_neurons_B

    ret_vals["rate_monitor_Z"] = rate_monitor_Z
    ret_vals["spike_monitor_Z"] = spike_monitor_Z
    ret_vals["voltage_monitor_Z"] = voltage_monitor_Z
    ret_vals["idx_monitored_neurons_Z"] = idx_monitored_neurons_Z

    ret_vals["rate_monitor_inhib"] = rate_monitor_inhib
    ret_vals["spike_monitor_inhib"] = spike_monitor_inhib
    ret_vals["voltage_monitor_inhib"] = voltage_monitor_inhib
    ret_vals["idx_monitored_neurons_inhib"] = idx_monitored_neurons_inhib

    return ret_vals




from neurodynex.tools import plot_tools
import matplotlib.pyplot as plt

for neverused in range(5):
	results =sim_decision_making_network_full()
	tA = results["rate_monitor_A"].t
	tA = [1000*tAA for tAA in tA]
	rA = results["rate_monitor_A"].smooth_rate(window="flat", width = 50 * b2.ms)
	tB = results["rate_monitor_B"].t
	tB = [1000*tBB for tBB in tB]
	rB = results["rate_monitor_B"].smooth_rate(window="flat", width = 50 * b2.ms)
	plt.plot(tA, rA, tB, rB)
	plt.xlabel('Time(ms)')
	plt.ylabel('Rate(Hz')
	plt.legend('Left', 'Right')
plt.show()
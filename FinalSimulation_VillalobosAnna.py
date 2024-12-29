from brian2 import *

# ------------------------------------------------------------- Parameters ------------------------------------------------------------
g_Ca = 4*mS/cm**2
g_K = 8*mS/cm**2
g_L = 2*mS/cm**2
V_Ca = 120*mV
V_K = -80*mV
V_L = -60*mV
phi = 1/(15.0*ms)
beta_1 = -1.2*mV
beta_2 = 18*mV
beta_3 = 10*mV
beta_4 = 17.4*mV
beta_5 = -1.2*mV
beta_6 = 18*mV
beta_7 = 10*mV
beta_8 = 17.4*mV
C1 = 1*uF/cm**2
C2 = 1*uF/cm**2
C3 = 1*uF/cm**2
C4 = 1*uF/cm**2
C5 = 1*uF/cm**2
C6 = 1*uF/cm**2
Iapp_1 = 50*uA/cm**2
Iapp_2 = 50*uA/cm**2
Iapp_5 = 10*uA/cm**2
Iapp_6 = 10*uA/cm**2
epsilon = 0.01 * mS/cm**2 # This is gap-junction coupling on the soma 
g_gap_a = 0.0000*mS/cm**2
g_gap_d = 0.00*mS/cm**2
g_c_d = 0.01*mS/cm**2 # Coupling Soma-Dendrite
g_c_a = g_c_d # Coupling Soma-Axon rn set to the same as dendrite
p_d = 0.1 # geometry parameter (relative size of soma vs dendrite)
p_a = p_d # geometry parameter (relative size of soma vs axon)

# ------------------------------------------------------------- Morphology ------------------------------------------------------------

# Same morphology for both neurons
#morpho = Soma(30*um)
#morpho.axon = Cylinder(diameter=1*um, length=300*um, n=1)
#morpho.dendrite = Cylinder(diameter=1*um, length=100*um, n=1)

gL = 0.002*siemens/cm**2
EL = -60*mV

# ---------------------------------------------------------------- Model ---------------------------------------------------------------

eqs = '''
dV1/dt = (-g_Ca*minf1*(V1-V_Ca)-g_K*n1*(V1-V_K)-g_L*(V1-V_L)+Iapp_1+epsilon*(V2-V1)+g_c_d*(V3-V1)/(1-p_d)+g_c_a*(V5-V1)/(1-p_a))/C1 : volt
dV2/dt = (-g_Ca*minf2*(V2-V_Ca)-g_K*n2*(V2-V_K)-g_L*(V2-V_L)+Iapp_2+epsilon*(V1-V2)+g_c_d*(V4-V2)/(1-p_d)+g_c_a*(V6-V2)/(1-p_a))/C2 : volt
dV3/dt = (-g_L*(V3 - V_L) + g_gap_d*(V4 - V3)+g_c_d*(V1-V3)/p_d)/C3 : volt
dV4/dt = (-g_L*(V4 - V_L) + g_gap_d*(V3 - V4)+g_c_d*(V2-V4)/p_d)/C4 : volt
dV5/dt = (-g_Ca*minf5*(V5-V_Ca)-g_K*n5*(V5-V_K)-g_L*(V5-V_L)+Iapp_5+g_gap_a*(V6-V5)+g_c_a*(V1-V5)/(p_a))/C5 : volt
dV6/dt = (-g_Ca*minf6*(V6-V_Ca)-g_K*n6*(V6-V_K)-g_L*(V6-V_L)+Iapp_6+g_gap_a*(V5-V6)+g_c_a*(V2-V6)/(p_a))/C6 : volt
dn1/dt = (ninf1-n1)/taun1 : 1
dn2/dt = (ninf2-n2)/taun2 : 1
dn5/dt = (ninf5-n5)/taun5 : 1
dn6/dt = (ninf6-n6)/taun6 : 1
minf1 = 0.5*(1 + tanh((V1 - beta_1)/beta_2)) : 1
minf2 = 0.5*(1 + tanh((V2 - beta_1)/beta_2)) : 1
minf5 = 0.5*(1 + tanh((V5 - beta_1)/beta_2)) : 1
minf6 = 0.5*(1 + tanh((V6 - beta_1)/beta_2)) : 1
taun1 = (phi*cosh((V1 - beta_3)/(2*beta_4)))**-1 : second
taun2 = (phi*cosh((V2 - beta_3)/(2*beta_4)))**-1 : second
taun5 = (phi*cosh((V5 - beta_3)/(2*beta_4)))**-1 : second
taun6 = (phi*cosh((V6 - beta_3)/(2*beta_4)))**-1 : second
ninf1 = 0.5*(1 + tanh((V1 - beta_3)/beta_4)) : 1
ninf2 = 0.5*(1 + tanh((V2 - beta_3)/beta_4)) : 1
ninf5 = 0.5*(1 + tanh((V5 - beta_3)/beta_4)) : 1
ninf6 = 0.5*(1 + tanh((V6 - beta_3)/beta_4)) : 1
'''

# ---------------------------------------------------------- Creation neurons ----------------------------------------------------------

# Neuron Group
N = 1
neurons = NeuronGroup(N, eqs, method='exponential_euler')

# Somas
neurons.V1 = 40*mV
neurons.n1 = 0.04
neurons.V2 = -40*mV
neurons.n2 = 0.074

# Dendrites
neurons.V3 = -60*mV
neurons.V4 = -60*mV

# Axons
neurons.V5 = 60*mV
neurons.n5 = 0.06
neurons.V6 = -60*mV
neurons.n6 = 0.06


# ---------------------------------------------------------- Monitors ----------------------------------------------------------
V1_monitor = StateMonitor(neurons, 'V1', record=True)
V2_monitor = StateMonitor(neurons, 'V2', record=True)
V3_monitor = StateMonitor(neurons, 'V3', record=True)
V4_monitor = StateMonitor(neurons, 'V4', record=True)
V5_monitor = StateMonitor(neurons, 'V5', record=True)
V6_monitor = StateMonitor(neurons, 'V6', record=True)
n1_monitor = StateMonitor(neurons, 'n1', record=True)
n2_monitor = StateMonitor(neurons, 'n2', record=True)

# ---------------------------------------------------------- Simulation ----------------------------------------------------------

run(5000 * ms)

print("SIM DONE")

# ---------------------------------------------------------- Plotting ----------------------------------------------------------
fig, axs = plt.subplots(3, 2, figsize=(15, 10))
axs[0, 0].plot(V1_monitor.t/ms, V1_monitor.V1[0]/mV,label='V1')
axs[0, 0].plot(V2_monitor.t/ms, V2_monitor.V2[0]/mV,  label='V2')
axs[0, 0].set(xlabel='Time (ms)', ylabel='Voltage (mV)', title='Soma Voltages')
axs[0, 0].legend()

axs[0, 1].plot(V3_monitor.t/ms, V3_monitor.V3[0]/mV, label='V3')
axs[0, 1].plot(V4_monitor.t/ms, V4_monitor.V4[0]/mV, label='V4')
axs[0, 1].set(ylim=(-70,-50), xlabel='Time (ms)', ylabel='Voltage (mV)', title='Dendrite Voltages')
axs[0, 1].legend()

axs[1, 0].plot(V1_monitor.t/ms, V1_monitor.V1[0]/mV, label='V1')
axs[1, 0].plot(V2_monitor.t/ms, V2_monitor.V2[0]/mV, label='V2')
axs[1, 0].set(xlabel='Time (ms)', ylabel='Voltage (mV)', title='Selected Interval Voltages')
axs[1, 0].set_xlim(100, 200)
axs[1, 0].legend()

axs[1, 1].plot(V5_monitor.t/ms, V5_monitor.V5[0]/mV, label='V5')
axs[1, 1].plot(V6_monitor.t/ms, V6_monitor.V6[0]/mV, label='V6')
axs[1, 1].set(xlabel='Time (ms)', ylabel='Voltage (mV)', title='Axon Voltages')
axs[1, 1].legend()

composed_signal = V1_monitor.V1[0] + V2_monitor.V2[0]
axs[2, 0].plot(V1_monitor.t/ms, composed_signal/mV, label='Composed Signal')
axs[2, 0].set(xlabel='Time (ms)', ylabel='Voltage (mV)', title='Composed Signal')
axs[2, 0].legend()

# Compute periodograms
from scipy.signal import periodogram

def compute_periodogram(signal, dt):
    f, Pxx = periodogram(signal, 1/dt)
    return f, Pxx

f1, Pxx1 = compute_periodogram(V1_monitor.V1[0], defaultclock.dt)
f2, Pxx2 = compute_periodogram(V2_monitor.V2[0], defaultclock.dt)

axs[2, 1].plot(f1, Pxx1, label='V1 Periodogram')
axs[2, 1].plot(f2, Pxx2, label='V2 Periodogram')
axs[2, 1].set(xlabel='Frequency (Hz)', ylabel='Power', title='Periodograms')
axs[2, 1].legend()

plt.tight_layout()
plt.show()

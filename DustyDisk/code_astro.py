import numpy as np
import matplotlib.pyplot as plt

# everying in cgs

G = 6.6743e-8         #[cm^3 g^-1 s^-2]
M_star = 1.99e33      
rho_s = 3.0           # dust density (material) [g/cm^3]

AU = 1.5e13
r_min = 0.1 * AU     #0.1 AU
r_max = 100 * AU     #100 AU
N = 1000
r = np.linspace(r_min, r_max, N)
dr = np.gradient(r)

#can add pressure bump at 30 AU ?
r0 = 30 *AU
sigma0 = 200     # [g/cm^2] at 1 AU or smth
Sigma_g = (sigma0 * (r / AU)**(-1.0)) * (1 + 1 * np.exp(-((r - r0) / (5*AU))**2))
           # normal power law            # bump

# temp and cs
T0 = 280     # K at 1 AU
mu = 2.34
k_B = 1.3807e-16
m_H = 1.67e-24
cs = np.sqrt(k_B * T0 * (r / (AU))**(-0.5) / (mu * m_H))  # [cm/s]

#kep omega and vk
Omega_K = np.sqrt(G * M_star / r**3)
v_K = Omega_K * r

# pressure
P = Sigma_g * cs * Omega_K / np.sqrt(2 * np.pi)   # dyn/cm^2

# eta  , such that v gas = vk (1 - eta) for small eta (taylor expansion)
dlnP = np.gradient(np.log(P), np.log(r))
eta = -0.5 * (cs / v_K)**2 * dlnP


sizes = [1e-3, 1e-2, 0.5, 2.0,10.0]  #the sizes are making sense?
colors = ['b', 'g', 'orange', 'r','k']
labels = [f'{s:.0e} cm' for s in sizes]



t_final = 1e12  # approx 30,000 yrs   (more? less?)
dt = 1e9        # approx 30 yrs   
Nt = int(t_final // dt)


dust_densities = []
drift_velocities = []

for s in sizes:
    
    St = (np.pi / 2) * rho_s * s / Sigma_g

    
    v_drift = -2 * eta * v_K * St / (1 + St**2)

    
    Sigma_d = np.ones_like(r)

    for _ in range(Nt):
        F = Sigma_d * v_drift
        dSigma_dt = -1 / r * np.gradient(r * F, r)
        Sigma_d += dSigma_dt * dt
        Sigma_d = np.maximum(Sigma_d, 1e-20)  # avoid negative values

    Sigma_d /= np.max(Sigma_d)

    dust_densities.append(Sigma_d)
    drift_velocities.append(v_drift)

fig, axs = plt.subplots(4, 1, figsize=(10, 15), sharex=True)

# gas properties
ax = axs[0]

ax.plot(r / AU, Sigma_g, lw=2, label='Gas Density')
ax.set_yscale('log')

ax.tick_params(axis='y')
ax.plot(r / AU, P,  lw=2, ls='--', label='Gas Pressure')


ax.grid(which='both', linestyle=':', linewidth=0.7)
ax.legend(loc='lower left', fontsize=10)


# dust density and stokes parameters

ax = axs[1]
for Sigma_d, color, label in zip(dust_densities, colors, labels):
    ax.loglog(r /AU, Sigma_d, color=color, lw=2, label=label)
ax.set_ylabel('norm dust den', fontsize=12)
ax.set_yscale('linear')
ax.legend(title="Grain size", fontsize=10, title_fontsize=11)
ax.grid(True, linestyle=':', linewidth=0.7)

#Drift velocity for different grain sizes
ax = axs[2]
for v_drift, color, label in zip(drift_velocities, colors, labels):
    ax.plot(r / AU, v_drift, color=color, lw=2, label=label)
ax.set_xlabel("Radius [AU]", fontsize=12)
ax.set_ylabel("Drift Velocity [cm/s]", fontsize=12)
ax.legend(title="Grain size", fontsize=10, title_fontsize=11)
ax.grid(True, linestyle=':', linewidth=0.7)


ax = axs[3]
ax.plot(r/AU , eta , label='eta')
ax.legend()
ax.grid(True, linestyle=':', linewidth=0.7)


for ax in axs:
    ax.set_xscale('log')
    ax.set_xlim(r_min / AU, r_max / AU)
    ax.minorticks_on()

plt.tight_layout()
plt.show()

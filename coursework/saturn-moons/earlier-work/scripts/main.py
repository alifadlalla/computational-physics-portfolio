from input_handler import load_input
from modules import Body, Simulation, NumericalMethod, Acceleration
import numpy as np

# Load input data
#input_file = 'input_data.json'
#bodies, sim_params = load_input(input_file)


# Unpack simulation parameters
#nbodies = sim_params['nbodies']
#dt = sim_params['dt']
#nt = sim_params['nt']


# Define celestial bodies in the Saturnian system
Saturn = Body('Saturn',
              37931206.234,  # GM (km^3/s^2),
              58232, # mean radius (km)
              np.array([0., 0., 0.]),  # Position (km)
              np.array([0., 0., 0.]),  # Velocity (km/s)
              'red')

Mimas = Body('Mimas',
             2.503489,
             198.8,
             np.array([-1.781644589918471E+05, 4.978604034734525E+04, 5.015471754576375E+03]),
             np.array([-3.101983306920243E+05, -1.200735430939711E+06, 4.624084160095980E+03])/(24*60*60),
             'orange')

Tethys = Body('Tethys',
              41.21,
              536.3,
              np.array([2.874974592673592E+05, 6.446905997181479E+04, 5.205343659123462E+03]),
              np.array([-2.142954962270242E+05, 9.570034498442297E+05, -7.058677519081212E+03])/(24*60*60),
              'green')

Titan = Body('Titan',
             8978.14,
             2575.5, # radius (Km)
             np.array([-7.931453936081687E+05, 9.251201308131195E+05, -7.317676996395952E+03]),
             np.array([-3.573469027363011E+05, -3.245883957332915E+05, -1.661041185020754E+03])/(24*60*60),
             'blue')

saturnian_system = [Saturn, Mimas, Tethys, Titan]
nbodies = 4
#dt = 2000
nt = 5000
dt = (24.8*24*60*60)/nt
print()
# Run simulations
#rk4_sim = Simulation(saturnian_system, dt, nt, method=NumericalMethod.runge_kutta, perturb=Acceleration.perturb, oblat=Acceleration.oblat)
rk4_sim = Simulation(saturnian_system, dt, nt, method=NumericalMethod.runge_kutta, perturb=Acceleration.perturb, oblat=Acceleration.oblat)
lf_sim = Simulation(saturnian_system, dt, nt, method=NumericalMethod.leapfrog,perturb=Acceleration.perturb, oblat=Acceleration.oblat)

#print(rk4_sim)
# Visualization or analysis here
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def plot_3d(sim, t):
    """3D plot of the simulation at time step `t`."""
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((1, 1, 1))

    c = [body.color for body in sim.bodies]
    b = [body.name for body in sim.bodies]
    s = [body.radius for body in sim.bodies]
    s = s * 100 / np.linalg.norm(s)

    x, y, z = np.rollaxis(sim.pos, 2)
    
    # Plot central body
    ax.scatter(x[t, 0], y[t, 0], z[t, 0], c=c[0], s=s[0]*50, label=b[0])

    # Plot paths and bodies
    for i in range(1, nbodies):
        ax.plot(x[:t, i], y[:t, i], z[:t, i], c=c[i], label=b[i])
        ax.scatter(x[t, i], y[t, i], z[t, i], c=c[i], s=s[i])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()
    plt.show()


def plot_ell(rk4_orb, frog_orb, syst):
    """Plot orbital elements."""
    fig, axs = plt.subplots(3, 6, sharex=True, figsize=(12, 6))
    lab = [r'$a(t)$', r'$e(t)$', r'$I(t)$', r'$☊(t)$', r'$\varpi(t)$', r'$\lambda(t)$']
    bcl = [body.color for body in syst][1:]
    bnm = [body.name for body in syst][1:]

    for i in range(nbodies-1):  # Iterate over bodies
        for j in range(6):  # Iterate over orbital elements
            axs[i, j].plot(rk4_orb[i, j], c=bcl[i])
            axs[i, j].plot(frog_orb[i, j], ls='--', c=bcl[i])
            axs[-1, j].set(xlabel=lab[j])
            axs[i, 0].set(ylabel=bnm[i])
            axs[i, j].yaxis.set_major_formatter(FormatStrFormatter('%.2e'))

    fig.suptitle('Orbital Elements (continuous: RK4, dashed: Leapfrog)')
    plt.tight_layout()
    plt.savefig('orbital_elements.png', dpi=300)
    plt.show()


def plot_libration_angles(orb_rk4):
    """Plot libration angles."""
    l1, l2, l3 = orb_rk4[(1, 2, 3), 5]
    w1, w2 = orb_rk4[(1, 2), 4]

    phi1 = l1 - 2 * l2 + w1
    phi2 = l1 - 2 * l2 + w2
    phi3 = l2 - 2 * l3 + w2

    phis = np.zeros_like(phi2 - phi3)
    for i, phi in enumerate(phi2 - phi3):
        while phi > np.pi:
            phi -= 2 * np.pi
        while phi < -np.pi:
            phi += 2 * np.pi
        phis[i] = phi

    plt.plot(phis, label=r'$\phi_3$')
    plt.ylim([-np.pi, np.pi])
    plt.legend()
    plt.show()


# Example simulation parameters
#if __name__ == "__main__":
    

# Generate plots
plot_3d(rk4_sim, 500)  # Adjust time step as needed
# You would need to generate orbital elements for plot_ell
plot_ell(rk4_sim.get_orbital_elements(), lf_sim.get_orbital_elements(), saturnian_system)
#print(rk4_sim.get_orbital_elements())
#plot_libration_angles(rk4_sim.get_orbital_elements())


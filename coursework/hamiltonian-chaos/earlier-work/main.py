import numpy as np
import matplotlib.pyplot as plt

def chirikov_map(I, theta, K, num_iterations):
    I_values, theta_values = [], []

    for _ in range(num_iterations):
        I_values.append(I)
        theta_values.append(theta)

        I_next = (I + K * np.sin(theta)) % (2 * np.pi)
        theta_next = (theta + I_next) % (2 * np.pi)

        I, theta = I_next, theta_next

    return I_values, theta_values

# Set parameters
K_values = [0.1, 0.25, 0.5, 0.9, 2, 4]  # Different values of K
#K_values = [0.1]  # Different values of K
M = 20  # Size of initial grid
N = 2000  # Number of iterations
L = 1000  # Size of array for phase space


scales=[1.0]
#scales = [1.0, 0.5, 0.25, 0.125, 0.001]  # Different scales

# Create initial grid excluding 2pi
I_initial = np.linspace(0, 2*np.pi-0.001, M)
theta_initial = np.linspace(0, 2*np.pi-0.001, M)

# Create array for phase space
phase_space = np.zeros((L, L))

# Iterate over different values of K
for K in K_values:
    # Iterate over different scales
    for scale in scales:
        I_initial = np.linspace(0, 2*np.pi*scale, M)
        theta_initial = np.linspace(0, 2*np.pi*scale, M)
        # Iterate over initial grid
        for I_init in I_initial:
            for theta_init in theta_initial:
                # Generate the map
                I_values, theta_values = chirikov_map(I_init, theta_init, K, N)

                # Convert to indices in the phase space array
                I_indices = (np.array(I_values) / (2*np.pi) * L).astype(int) % L
                theta_indices = (np.array(theta_values) / (2*np.pi) * L).astype(int) % L

                # Light up the visited cells
                phase_space[I_indices, theta_indices] = 1

        # Plot the phase space
        #plt.figure(figsize=(6, 6), dpi=200)
        fig = plt.figure(figsize=(6, 6), dpi=3000)
        plt.imshow(phase_space, extent=[0, 2*np.pi*scale-0.0001, 0, 2*np.pi*scale-0.0001], origin='lower')
        
        # Customize x and y ticks in terms of pi
        plt.xticks([0, np.pi, 2*np.pi], ['0', '$\pi$', '$2\pi$'], fontsize=5)
        plt.yticks([0, np.pi, 2*np.pi], ['0', '$\pi$', '$2\pi$'], fontsize=5)
        
        plt.title('Phase Space for K = {}'.format(K), fontsize=10)
        plt.xlabel('I', fontsize=10)
        plt.ylabel(r'$\theta$', fontsize=10)
        k_str = str(K).replace('.', '_')
        plt.savefig(f'K_{k_str}.pdf', format='pdf', dpi=3000, bbox_inches='tight')
        print('Done')
        #plt.show()

        # Reset the phase space for the next value of K and scale
        phase_space.fill(0)


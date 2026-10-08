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

def calculate_lyapunov_exponent(I, theta, K, num_iterations):
    lyapunov_sum = 0

    for _ in range(num_iterations):
        I_next, theta_next = chirikov_map(I, theta, K, 1)
        
        # Calculate the derivative of the map
        df_dI = 1 + K * np.cos(theta)

        # Update the Lyapunov sum
        lyapunov_sum += np.log(np.abs(df_dI))

        I, theta = I_next[-1], theta_next[-1]

    return lyapunov_sum / num_iterations

def iterate(M,N,L,K,ext='svg', lyapov=False):

    # Create initial grid excluding 2pi
    I_initial = np.linspace(0, 2*np.pi-0.001, M)
    theta_initial = np.linspace(0, 2*np.pi-0.001, M)

    # Create array for phase space
    phase_space = np.zeros((L, L))

    # Iterate over initial grid
    for I_init in I_initial:
        for theta_init in theta_initial:
            # Generate the map
            I_values, theta_values = chirikov_map(I_init, theta_init, K, N)

            # Convert to indices in the phase space array
            I_indices = (np.array(I_values) / (2*np.pi) * L).astype(int) % L
            theta_indices = (np.array(theta_values) / (2*np.pi) * L).astype(int) % L

            

            if lyapov:
                # Calculate Lyapunov exponent for the cell
                lyapunov_exponent = calculate_lyapunov_exponent(I_init, theta_init, K, N)
                # Light up the visited cells with color based on Lyapunov exponent
                phase_space[I_indices, theta_indices] = lyapunov_exponent
            else:
                # Light up the visited cells
                phase_space[I_indices, theta_indices] = 1


    # Plot the phase space
    DPI = 300
    fig = plt.figure(figsize=(6, 6), dpi=DPI)
    # Customize x and y ticks in terms of pi
    plt.xticks([0, np.pi, 2*np.pi], ['0', '$\pi$', '$2\pi$'], fontsize=5)
    plt.yticks([0, np.pi, 2*np.pi], ['0', '$\pi$', '$2\pi$'], fontsize=5)
    
    plt.title('Phase Space for K = {}'.format(K), fontsize=10)
    plt.xlabel('I', fontsize=10)
    plt.ylabel(r'$\theta$', fontsize=10)
    k_str = str(K).replace('.', '_')
    if lyapov:
        plt.imshow(phase_space, extent=[0, 2*np.pi, 0, 2*np.pi], origin='lower', cmap='coolwarm')
        plt.colorbar(label='Mean Lyapunov Exponent', shrink=0.75)
        if ext=='pdf':
            plt.savefig(f'K_{k_str}_lyapov.pdf', format='pdf', dpi=DPI, bbox_inches='tight')
        elif ext =='svg':
            plt.savefig(f'K_{k_str}_lyapov.svg', format='svg', dpi=DPI, bbox_inches='tight')
    else:
        plt.imshow(phase_space, extent=[0, 2*np.pi-0.0001, 0, 2*np.pi-0.0001], origin='lower')
        if ext=='pdf':
            plt.savefig(f'K_{k_str}.pdf', format='pdf', dpi=DPI, bbox_inches='tight')
        elif ext =='svg':
            plt.savefig(f'K_{k_str}.svg', format='svg', dpi=DPI, bbox_inches='tight')
    plt.show()



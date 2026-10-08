import numpy as np
import matplotlib.pyplot as plt

def plot_speedup(num_threads, sequential_time, parallel_times, label):
    """Plots the speedup curve."""
    speedup = sequential_time / parallel_times
    plt.plot(num_threads, speedup, marker='o', linestyle='-', label=label)

def plot_amdahl(num_threads, p):
    """Plots Amdahl's Law speedup."""
    speedup = 1 / ((1 - p) + p / num_threads)
    plt.plot(num_threads, speedup, linestyle='--', label=f"Amdahl's Law (P={p:.3f})")

if __name__ == "__main__":
    
    num_threads = np.array([1, 2, 4, 6, 8, 16, 24, 25, 32, 40, 48])
    
    # Data for the first version
    sequential_time1 = 74228
    parallel_times1 = np.array([143713, 74309, 39066, 27404, 21879, 12459, 9269, 14222, 11890, 9430, 10451])

    # Data for the second version
    sequential_time2 = 86971
    parallel_times2 = np.array([145478, 74628, 39541, 28425, 22049, 12867, 10970, 14955, 13162, 13200, 11674])

    # Calculate parallel fraction (P) using Amdahl's Law (example, adjust if needed)
    # Using the max speedup of the first version with 24 threads:
    speedup_24_1 = sequential_time1 / parallel_times1[6] # index 6 correspond to 24 threads
    p1 = (1 - (1/speedup_24_1)) / (1 - (1/24))
    #print(f"Parallel fraction for the first version: {p1:.3f}")

    # Using the max speedup of the second version with 24 threads:
    speedup_24_2 = sequential_time2 / parallel_times2[6] # index 6 correspond to 24 threads
    p2 = (1 - (1/speedup_24_2)) / (1 - (1/24))
    #print(f"Parallel fraction for the second version: {p2:.3f}")

    # Plot the speedup curves
    plt.figure(figsize=(10, 6))

    plot_speedup(num_threads, sequential_time1, parallel_times1, label="Parallelized Accel.")
    plot_speedup(num_threads, sequential_time2, parallel_times2, label="Parallelized Accel., Pos., Vel., Ec.")
    plot_amdahl(num_threads, p1)
    #plot_amdahl(num_threads, p2) same as p1


    plt.xlabel("Number of Threads")
    plt.ylabel("Speedup")
    plt.title("Speedup Comparison")
    plt.xticks(np.unique(np.concatenate((num_threads, num_threads))))
    plt.grid(True)
    plt.legend()
    #plt.savefig("speedup_comparison.pdf")
    plt.show()
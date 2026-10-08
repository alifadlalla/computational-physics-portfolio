import yaml
import numpy as np
import matplotlib.pyplot as plt

def read_input_yaml(filename="../input.yaml"):
    """Reads parameters from a YAML input file."""
    try:
        with open(filename, 'r') as f:
            data = yaml.safe_load(f)
            return data
    except FileNotFoundError:
        print(f"Error: Input file '{filename}' not found.")
        return None
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")
        return None

def load_energy_data(file_name,nstep):
    """Load energy data from a binary file."""
    try:
        with open(file_name, 'rb') as f:
            data = np.fromfile(f, dtype=np.float64).reshape(nstep, 2)
            return data
    except IOError:
        print(f"Error: Could not read file {file_name}.")
        return None
    except ValueError as e:
        print(f"Error reading energy file: {e}")
        return None

def plot_energies(data):
    """Plot kinetic and potential energies over time."""
    if data is None:
        print("No data to plot.")
        return

    potential_energy = data[:, 0]
    kinetic_energy = data[:, 1]
    time_steps = np.arange(len(potential_energy))
    total_energy = potential_energy + kinetic_energy

    plt.figure(figsize=(10, 6))
    plt.plot(time_steps, potential_energy, label="Potential Energy", color="blue", lw=2)
    plt.plot(time_steps, kinetic_energy, label="Kinetic Energy", color="red", lw=2)
    plt.plot(time_steps, total_energy, label="Total Energy", color="green", lw=2, linestyle="--")
    print("percentage=", abs((max(total_energy) - min(total_energy))*100/(max(total_energy))), "%")
    plt.title("Energy Evolution Over Time", fontsize=16)
    plt.xlabel("Time Step", fontsize=14)
    plt.ylabel("Energy", fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.savefig("energy.pdf")
    plt.show()


if __name__ == "__main__":
    input_data = read_input_yaml()
    if input_data:
        nstep = input_data.get('nstep')
    else:
        exit()

    energy_file = "energy.bin"
    energy_data = load_energy_data(energy_file,nstep)

    if energy_data is not None :
        plot_energies(energy_data)
    else:
        print("Could not load data. Check if the files exist and are correctly formatted.")
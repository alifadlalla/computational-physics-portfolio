import yaml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # Import for 3D plotting

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


def load_position_data(file_name,npart,nstep):
    try:
        with open(file_name, 'rb') as f:
            positions = np.fromfile(f, dtype=np.float64).reshape(nstep+1, 3, npart)
            return positions
    except IOError:
        print(f"Error: Could not read file {file_name}.")
        return None
    except ValueError as e:
        print(f"Error reading position file: {e}")
        return None

def animate_particles(positions, npart, nstep, figs=(800, 600)):
    """Animates particle movement in 3D with improved visuals and sphere outline.

    Args:
        positions (numpy.ndarray): Array of particle positions (nstep+1, 3, npart).
        npart (int): Number of particles.
        nstep (int): Number of time steps.
        figs (tuple): Figure size (width, height) in pixels.
    """

    dpi_v = 100
    sfig = (figs[0] / dpi_v, figs[1] / dpi_v)

    fig = plt.figure(figsize=sfig, facecolor='k', num=1)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.axes.set_xlim3d(left=-1.5, right=1.5)
    ax.axes.set_ylim3d(bottom=-1.5, top=1.5)
    ax.axes.set_zlim3d(bottom=-1.5, top=1.5)
    ax.tick_params(axis='x', colors='w', which='both')
    ax.tick_params(axis='y', colors='w', which='both')
    ax.tick_params(axis='z', colors='w', which='both')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] = (0.25, 0.25, 0.25, 1)
    ax.yaxis._axinfo["grid"]['color'] = (0.25, 0.25, 0.25, 1)
    ax.zaxis._axinfo["grid"]['color'] = (0.25, 0.25, 0.25, 1)
    ax.set_facecolor("k")
    plt.tight_layout(pad=1)

    # Draw representative circles of the sphere
    theta = np.linspace(0, 2 * np.pi, 100)

    # xy-plane circle
    x = np.cos(theta)
    y = np.sin(theta)
    z = np.zeros_like(theta)
    ax.plot(x, y, z, color="r",linestyle="--", alpha=0.5)

    # xz-plane circle
    x = np.cos(theta)
    z = np.sin(theta)
    y = np.zeros_like(theta)
    ax.plot(x, y, z, color="r",linestyle="--", alpha=0.5)

    # yz-plane circle
    y = np.cos(theta)
    z = np.sin(theta)
    x = np.zeros_like(theta)
    ax.plot(x, y, z, color="r",linestyle="--", alpha=0.5)


    graph = ax.scatter(positions[0, 0, :], positions[0, 1, :], positions[0, 2, :], s=2, c='w', linewidths=0)

    ax.view_init(elev=20, azim=20)  # Set the initial view

    def update_graph(t):
        graph._offsets3d = (positions[t, 0, :], positions[t, 1, :], positions[t, 2, :])
        return graph,
    
    ani = animation.FuncAnimation(fig, update_graph, frames=nstep + 1, interval=1, blit=False)
    plt.show()

def plot_initial_positions(positions):
    """Plots only the initial particle positions.

    Args:
        positions (numpy.ndarray): The full position data (nstep+1, 3, npart).
    """
    if positions is None:
        print("No position data to plot.")
        return

    initial_positions = positions[0, :, :]  # Extract initial positions (time step 0)
    print(np.all(np.array(initial_positions) < 1))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(initial_positions[0, :], initial_positions[1, :], initial_positions[2, :], s=1, c='b')

    # Draw a wireframe sphere for reference
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(x, y, z, color="r", alpha=0.3)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Initial Particle Positions')
    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_zlim([-1.2, 1.2])
    plt.show()

if __name__ == "__main__":
    input_data = read_input_yaml()
    if input_data:
        npart = input_data.get('npart')
        nstep = input_data.get('nstep')
    else:
        exit()

    position_file = "position.bin"
    position_data = load_position_data(position_file,npart,nstep)

    if position_data is not None:
        animate_particles(position_data, npart,nstep)
        #plot_initial_positions(position_data) # to confirm that the particles are in the sphere of r=1
    else:
        print("Could not load data. Check if the files exist and are correctly formatted.")
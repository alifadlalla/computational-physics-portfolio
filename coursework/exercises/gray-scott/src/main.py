
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Assuming 'data_u.txt' and 'data_v.txt' are the files with values for u and v respectively
u_data = np.loadtxt('output_u.txt')
v_data = np.loadtxt('output_v.txt')

# Get the number of timesteps (lines in the file)
timesteps = u_data.shape[0]

# Reshape each line into a grid (adjust grid size as needed, e.g., 64x64)
grid_size = (1024, 1024)  # Example grid size, modify according to your data
u_data = u_data.reshape((timesteps, *grid_size))
v_data = v_data.reshape((timesteps, *grid_size))

# Set up the figure and axis
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Initialize the plots
u_plot = ax[0].imshow(u_data[0], cmap='viridis', vmin=u_data.min(), vmax=u_data.max())
v_plot = ax[1].imshow(v_data[0], cmap='plasma', vmin=v_data.min(), vmax=v_data.max())

# Set titles
ax[0].set_title('U Concentration')
ax[1].set_title('V Concentration')

# Function to update the plot at each frame
def update(frame):
    u_plot.set_data(u_data[frame])
    v_plot.set_data(v_data[frame])
    return u_plot, v_plot

# Create the animation
# Create the animation
anim = FuncAnimation(fig, update, frames=timesteps, interval=100, blit=True)

# Show the animation
plt.show()


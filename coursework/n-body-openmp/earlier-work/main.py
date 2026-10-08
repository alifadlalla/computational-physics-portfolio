import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
print(os.getcwd())

# Step 1: Load the position and velocity data
positions = np.loadtxt('2D_evolution_r.dat')  # Each row: [x1, y1, x2, y2, ...]
velocities = np.loadtxt('2D_evolution_v.dat') # Each row: [vx1, vy1, vx2, vy2, ...]

# Step 2: Determine the number of particles and time steps
num_time_steps, num_columns = positions.shape
num_particles = num_columns // 2  # Each particle has x and y coordinates

# Reshape the data for easier indexing
positions = positions.reshape(num_time_steps, num_particles, 2)  # Shape: (time_steps, particles, 2)
velocities = velocities.reshape(num_time_steps, num_particles, 2)  # Shape: (time_steps, particles, 2)

# Step 3: Set up the plot
fig, ax = plt.subplots()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_title('2D Particle Evolution')
ax.set_xlabel('X Position')
ax.set_ylabel('Y Position')

# Initialize scatter plot and quiver (arrows)
scat = ax.scatter([], [], s=50, c='blue', label='Particles')
quiver = ax.quiver([], [], [], [], angles='xy', scale_units='xy', scale=1, color='red', label='Velocity')

# Step 4: Update function for animation
def update(frame):
    # Update positions and velocities for this frame
    x = positions[frame, :, 0]
    y = positions[frame, :, 1]
    vx = velocities[frame, :, 0]
    vy = velocities[frame, :, 1]

    # Update the scatter plot positions
    scat.set_offsets(np.c_[x, y])

    # Update the quiver plot for velocities
    quiver.set_offsets(np.c_[x, y])
    quiver.set_UVC(vx, vy)

    return scat, quiver

# Step 5: Create the animation
anim = FuncAnimation(fig, update, frames=num_time_steps, interval=50, blit=True)

# Step 6: Display the animation
plt.legend()
plt.show()

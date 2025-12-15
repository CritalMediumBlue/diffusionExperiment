import numpy as np
import pandas as pd

# Simulation parameters
n_particles = 10000
n_steps = 2000
L = 10.0
D_left = 1.0
D_right = 4.0
dt = 0.01

# Derived step sizes (discrete jump model)
# step_size = sqrt(2 * D * dt)
step_left = np.sqrt(2 * D_left * dt)
step_right = np.sqrt(2 * D_right * dt)

# Initialize particles uniformly (or at 0)
positions = np.random.uniform(-L, L, n_particles)

# Run simulation
for _ in range(n_steps):
    # Determine step size for each particle based on current position
    # Mask for left and right
    mask_left = positions < 0
    mask_right = ~mask_left
    
    # Generate random moves (+1 or -1) * step_size
    moves = np.random.choice([-1, 1], size=n_particles)
    
    # Apply step sizes
    # Note: This is the "Ito" interpretation mechanism: step size depends on *current* position
    current_steps = np.zeros(n_particles)
    current_steps[mask_left] = step_left
    current_steps[mask_right] = step_right
    
    # Update positions
    positions += moves * current_steps
    
    # Reflecting boundaries at -L and L
    # Simple reflection: if x > L, x = L - (x-L) = 2L - x
    # If x < -L, x = -L + (-L - x) = -2L - x
    
    # Using np.where for vectorization
    positions = np.where(positions > L, 2*L - positions, positions)
    positions = np.where(positions < -L, -2*L - positions, positions)

# Analyze results
# Count particles in left vs right (ignoring transients, assuming mixed)
# Use a histogram
bins = np.linspace(-L, L, 50)
hist, bin_edges = np.histogram(positions, bins=bins, density=True)

# Calculate ratio of mean density left vs right
mean_density_left = np.mean(hist[bin_edges[:-1] < 0])
mean_density_right = np.mean(hist[bin_edges[:-1] >= 0])
ratio = mean_density_left / mean_density_right

print(f"{mean_density_left=}")
print(f"{mean_density_right=}")
print(f"{ratio=}")
print(f"{D_right/D_left=}")

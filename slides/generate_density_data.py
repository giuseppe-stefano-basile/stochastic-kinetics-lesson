"""Reproduce the trajectory/forward-generator comparison in the slide deck.

Run from this directory: python generate_density_data.py
Requires NumPy and SciPy. Produces three CSV files. The TeX deck embeds the
resulting coordinates, so a LaTeX installation needs no Python to build it.
"""
from pathlib import Path
import numpy as np
from scipy.linalg import expm

D, beta, barrier, dt = 1.0, 1.0, 2.0, 0.002
left, right, n_grid, n_paths, seed = -2.0, 2.0, 81, 18000, 274
x = np.linspace(left, right, n_grid)
h = x[1] - x[0]
F = barrier * (x*x - 1)**2
Q = np.zeros((n_grid, n_grid))
for i in range(n_grid - 1):
    delta = beta * (F[i + 1] - F[i])
    Q[i, i + 1] = D / h**2 * np.exp(-delta / 2)
    Q[i + 1, i] = D / h**2 * np.exp(delta / 2)
np.fill_diagonal(Q, -Q.sum(axis=1))
initial = np.zeros(n_grid)
initial[np.argmin(abs(x + 1))] = 1
rng = np.random.default_rng(seed)
position = np.full(n_paths, -1.0)
samples = {100: 0.2, 500: 1.0, 2000: 4.0}
width = right - left
edges = np.r_[x[0]-h/2, (x[:-1]+x[1:])/2, x[-1]+h/2]
for step in range(1, max(samples) + 1):
    force = -4 * barrier * position * (position**2 - 1)
    position += D*beta*force*dt + np.sqrt(2*D*dt)*rng.normal(size=n_paths)
    # Fold the endpoint at both walls; exact reflection is approached as dt -> 0.
    position = left + np.abs((position-left+width) % (2*width)-width)
    if step in samples:
        time = samples[step]
        histogram = np.histogram(position, bins=edges)[0] / (n_paths*h)
        probability = initial @ expm(Q*time)
        assert np.isclose(probability.sum(), 1.0, atol=1e-10)
        data = np.column_stack((x, probability/h, histogram))
        np.savetxt(Path(f'density_t{time:g}.csv'), data, delimiter=',',
                   header='x,grid_density,trajectory_density', comments='')
        print(f't={time:g}: mass={probability.sum():.12f}, '
              f'mean |density difference|={np.mean(np.abs(data[:,1]-data[:,2])):.4f}')

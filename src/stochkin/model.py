"""Transparent 1D equilibrium diffusion models in dimensionless units."""
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.linalg import solve


def potential(x, barrier=2.0):
    x = np.asarray(x)
    return barrier * (x*x - 1)**2


def force(x, barrier=2.0):
    x = np.asarray(x)
    return -4 * barrier * x * (x*x - 1)


def mfpt_quadrature(x0, a, b, potential_fn=potential, diffusion=1.0,
                    beta=1.0, ngrid=4001):
    """Reflecting a, absorbing b; diffusion may be positive scalar or callable."""
    if not a <= x0 < b or ngrid < 3:
        raise ValueError("Require a <= x0 < b and ngrid >= 3")
    x = np.unique(np.r_[np.linspace(a, b, ngrid), x0])
    F = np.asarray(potential_fn(x)) + np.zeros_like(x)
    F = F - np.min(F)  # constant shifts have no effect
    D = np.asarray(diffusion(x) if callable(diffusion) else diffusion) + np.zeros_like(x)
    if np.any(D <= 0) or not np.all(np.isfinite(D)):
        raise ValueError("diffusion must be finite and positive")
    inner = cumulative_trapezoid(np.exp(-beta*F), x, initial=0)
    y = np.exp(beta*F)*inner/D
    i = np.searchsorted(x, x0)
    return float(trapezoid(y[i:], x[i:]))


def grid_generator(x, potential_fn=potential, diffusion=1.0, beta=1.0):
    """Row-generator Q on uniform reflecting grid; detailed balance exact at nodes.

    dp/dt = p Q for row probabilities, while Q acts on column observables.
    """
    x = np.asarray(x, dtype=float)
    if len(x) < 3 or diffusion <= 0 or not np.allclose(np.diff(x), x[1]-x[0]) or x[1] <= x[0]:
        raise ValueError("Need >=3 increasing uniformly spaced nodes and positive D")
    h = x[1]-x[0]
    F = np.asarray(potential_fn(x)) + np.zeros_like(x)
    Q = np.zeros((len(x), len(x)))
    for i in range(len(x)-1):
        d = beta*(F[i+1]-F[i])
        Q[i, i+1] = diffusion/h**2*np.exp(-d/2)
        Q[i+1, i] = diffusion/h**2*np.exp(d/2)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def grid_mfpt(Q, target=-1):
    """Expected continuous-time hitting time of target node from every node."""
    Q = np.asarray(Q)
    n = len(Q)
    target = target % n
    other = np.delete(np.arange(n), target)
    result = np.zeros(n)
    result[other] = solve(-Q[np.ix_(other, other)], np.ones(n-1))
    return result


def simulate_fpt(n_traj=2000, dt=0.002, max_time=150.0, x0=-1.0,
                 a=-2.0, b=1.0, D=1.0, beta=1.0, barrier=2.0, seed=42):
    """Return observed time and hit flag; censored times equal max_time.

    Euler–Maruyama endpoint crossings; left reflection uses repeated folding.
    Time is rounded upward to a full step. Decrease dt to assess crossing bias.
    """
    if n_traj < 1 or dt <= 0 or max_time <= 0 or D <= 0 or not a <= x0 < b:
        raise ValueError("Invalid simulation parameters")
    nsteps = int(round(max_time/dt))
    if nsteps < 1 or not np.isclose(nsteps*dt, max_time):
        raise ValueError("max_time must be an integer multiple of dt")
    rng = np.random.default_rng(seed)
    x = np.full(n_traj, x0, dtype=float)
    times = np.full(n_traj, max_time)
    hit = np.zeros(n_traj, dtype=bool)
    for step in range(1, nsteps+1):
        ids = np.flatnonzero(~hit)
        if not len(ids): break
        proposal = x[ids] + D*beta*force(x[ids], barrier)*dt + np.sqrt(2*D*dt)*rng.normal(size=len(ids))
        # A single reflecting wall: reflection of overshoot across a.
        proposal = a + np.abs(proposal-a)
        newly = proposal >= b
        times[ids[newly]] = step*dt
        hit[ids[newly]] = True
        x[ids] = proposal
    return times, hit

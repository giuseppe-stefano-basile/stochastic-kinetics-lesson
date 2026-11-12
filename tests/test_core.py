import numpy as np
from scipy.linalg import expm
from stochkin import mfpt_quadrature, grid_generator, grid_mfpt, potential, simulate_fpt


def flat(x): return np.zeros_like(x)


def test_flat_exact_and_diffusion_scaling():
    a, x, b, D = -2., -1., 1., 2.
    exact = ((b-a)**2-(x-a)**2)/(2*D)
    assert np.isclose(mfpt_quadrature(x,a,b,flat,D), exact, rtol=1e-6)
    assert np.isclose(mfpt_quadrature(x,a,b,potential,2), mfpt_quadrature(x,a,b,potential,1)/2)


def test_generator_and_grid_mfpt():
    x = np.linspace(-2,1,201)
    q = grid_generator(x, potential)
    pi = np.exp(-potential(x)); pi /= pi.sum()
    assert np.allclose(q.sum(axis=1),0,atol=1e-10)
    assert np.allclose(pi[:,None]*q, (pi[:,None]*q).T)
    assert np.allclose(pi @ q, 0, atol=1e-10)
    T = grid_mfpt(q)
    assert T[-1] == 0 and np.all(T[:-1]>0)
    theory = mfpt_quadrature(-1,-2,1)
    assert abs(np.interp(-1,x,T)/theory-1) < 0.04


def test_censoring_and_reproducibility():
    t,h = simulate_fpt(n_traj=30,max_time=0.02,dt=0.01,seed=7)
    t2,h2 = simulate_fpt(n_traj=30,max_time=0.02,dt=0.01,seed=7)
    assert np.array_equal(t,t2) and np.array_equal(h,h2)
    assert np.all(t[~h] == .02) and np.any(~h)



def test_variable_diffusion_flat_exact():
    # D(x)=1+x on [0,1], reflecting 0, absorbing 1.
    # T(0)=integral_0^1 y/(1+y)dy = 1-log(2).
    estimate = mfpt_quadrature(0., 0., 1., flat, lambda x: 1+x, ngrid=10001)
    assert np.isclose(estimate, 1-np.log(2), rtol=1e-7)


def test_fpt_simulation_agrees_with_quadrature_with_sampling_margin():
    # Fixed, moderate resolution; this exercises the simulation and the event flags.
    # Statistical tolerance also allows for endpoint-detection bias at finite dt.
    t, hit = simulate_fpt(n_traj=3000, dt=0.002, max_time=120., seed=2026)
    assert hit.all()
    exact = mfpt_quadrature(-1., -2., 1.)
    se = t.std(ddof=1)/np.sqrt(len(t))
    assert abs(t.mean()-exact) < 4*se + 0.20

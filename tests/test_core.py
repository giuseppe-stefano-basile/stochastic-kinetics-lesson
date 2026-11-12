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

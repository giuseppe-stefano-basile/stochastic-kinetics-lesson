# Derivations: stochastic dynamics and first passage

## 1. Probability and a stochastic trajectory

A process X_t is a family of random variables. A trajectory is one realization. The transition density p(x,t|x0,0) represents an ensemble conditional on X_0=x0. For Brownian motion W, increments over disjoint intervals are independent Gaussian random variables with mean 0 and variance equal to interval length. Thus ΔW=√Δt ξ, ξ~N(0,1). Informally (ΔW)² contributes at order Δt; individual trajectories are almost surely nondifferentiable. White noise is a distributional derivative, not an ordinary function.

## 2. Itô generator and forward equation

Assume a sufficiently smooth one-dimensional Itô diffusion, dX=a(X)dt+√(2D(X))dW, D>0. For a smooth test function f, Itô's formula is

    df = [a f' + D f'']dt + √(2D) f' dW.

Taking expectations yields d E[f(X_t)]/dt = E[L f(X_t)], where the **backward generator** Lf=af'+Df''. Write the expectation as ∫f p dx, integrate by parts twice, and require boundary terms to vanish under the appropriate domains. Then

    ∂t p = L* p = -∂x(ap)+∂x²(Dp) = -∂x J,
    J = ap - ∂x(Dp).

The forward equation evolves densities; L acts on observables of the initial point. The notation L and L* varies across textbooks. Under suitable regularity and boundary assumptions the same transition kernel solves the forward and backward equations in its final and initial coordinates, respectively.

## 3. Equilibrium and variable diffusion

At temperature T let β=(kBT)^−1. With D constant, a=−DβF' and

    ∂t p = ∂x[D(βF'p+p')],   J=−D(βF'p+p').

Zero-current boundaries give p_eq=Z^−1 exp(−βF), provided Z<∞. A stationary distribution need not have zero current on a periodic domain; here we assume equilibrium detailed balance. D is absent from p_eq, but multiplying D by c>0 accelerates the generator by c and divides first-passage times by c.

For spatially varying D, insisting on this equilibrium and Itô convention gives a=D'−DβF', because J=ap−(Dp)' must vanish for p∝exp(−βF). Therefore

    dX = [D'−DβF']dt + √(2D)dW,
    ∂t p = ∂x{D e^(−βF) ∂x[e^(βF)p]}.

This drift depends on the stochastic convention; equivalent Stratonovich representations require their own converted drift. Merely replacing constant D by D(x) in the drift gives the wrong stationary law. Units: [D]=[x]²/time, [F]=energy, and [βF]=1. On a nonlinear coordinate change both F and D transform; a free-energy profile is a density relative to its chosen coordinate measure.

## 4. Reflecting and absorbing boundaries

Let a<x<b, with a reflecting and b absorbing. In the forward equation reflection means J(a,t)=0; absorption is p(b,t)=0 for a regular diffusion. With an absorbing endpoint, p is a sub-probability density: its integral S(t)=P(τ_b>t) decreases according to dS/dt=−J(b,t). The nonnegative first-passage density is −S'(t)=J(b,t). The backward boundary conditions are T'(a)=0 and T(b)=0. Endpoint reflection used in the code approximates the continuous reflected process; finite-step endpoint detection misses some crossings.

## 5. Backward first-passage equation and integral

Define τ_b=inf{t≥0:X_t≥b} and T(x)=E_x[τ_b], assumed finite. Dynamic programming or Dynkin's formula gives L T=−1 with the boundary conditions above. For equilibrium diffusion,

    L T = e^(βF) d/dx [D e^(−βF) T'] = −1.

Multiplying by exp(−βF) and integrating from a to y, with T'(a)=0,

    D(y)e^(−βF(y))T'(y) = − ∫_a^y exp(−βF(z)) dz.

Integrating T' from x to b and using T(b)=0 yields

    T(x)= ∫_x^b dy [e^(βF(y))/D(y)] ∫_a^y dz e^(−βF(z)).

For F=0 and D constant this reduces to ((b−a)²−(x−a)²)/(2D). For a random nonnegative τ, Eτ=∫_0^∞P(τ>t)dt, including infinite values. Its empirical finite-time integral is a restricted mean, not necessarily the full MFPT. A fixed initial point, an equilibrium distribution within a basin, and a quasi-stationary starting distribution give different averages.

## 6. Committor and reaction rates

For absorbing sets at a and b, q(x)=P_x(τ_b<τ_a) satisfies Lq=0, q(a)=0, q(b)=1. With variable D,

    q(x)= [∫_a^x e^(βF(y))/D(y) dy] / [∫_a^b e^(βF(y))/D(y) dy].

1/T is the reciprocal mean for the specified first-passage experiment; it is not automatically an equilibrium two-state rate. A simple rate interpretation needs appropriate state definitions and suitable timescale separation, and depends on the initial distribution and treatment of recrossings. First passage to a barrier top and arrival in a product basin are distinct events.

## 7. Grid generator and MSM connection

For uniform nodes x_i, spacing h, define Q_{i,i+1}=D/h² exp[−β(F_{i+1}−F_i)/2] and the reverse rate with opposite sign. Set Q_ii=−sum_{j≠i} Q_ij, omitting jumps across reflecting endpoints. Rows sum to zero, off-diagonals are nonnegative, and π_i∝exp(−βF_i) obeys π_i Q_ij=π_j Q_ji. For row probability vectors p, dp/dt=pQ and P(t)=exp(tQ); for an observable column f, df/dt=Qf. Removing an absorbing target's row and column gives −Q_SS T=1. Finite-grid reflecting endpoints approximate the continuum boundary and require a grid convergence study. A fitted MSM transition matrix at finite lag is a separate approximation, not automatically this generator.

## 8. Molecular interpretation and limitations

For collective variable s(R), p_eq(s)=∫δ(s−s(R))ρ_eq(R)dR and F(s)=−kBT ln p_eq(s)+C. Integrating out coordinates generally introduces memory and coordinate-dependent transport. A 1D Markov diffusion is a testable approximation, not a consequence of knowing F(s). Check lag-time dependence, Markov/Chapman–Kolmogorov tests, conditional relaxation, first-passage predictions, and sensitivity to state boundaries; distinguish force-field error, uncertainty in F and D, and model error. Recovering a free energy from metadynamics does not automatically recover unbiased trajectory kinetics.

# Instructor guide

## Preparation

Install the package and run both notebooks from start to finish. Students need elementary calculus (integration by parts and Taylor expansions), basic probability, arrays and plotting. The course does not require prior stochastic calculus: introduce Wiener increments and Itô's formula before using them. Explain that a proof of existence/uniqueness of the SDE is beyond scope; the manipulations use smooth coefficients, suitable boundary regularity and integrability.

## Predictions and expected answers

1. Brownian motion: mean 0, variance 2Dt. The variance, not the standard deviation, is linear in t. The finite-ensemble curves fluctuate.
2. Double well: equilibrium density is proportional to exp(−βF). Increasing D at fixed F changes temporal sampling, not equilibrium; finite simulations can be misleading.
3. FPT: the target is an event definition; survival decreases, the finite-horizon restricted mean underestimates the full MFPT when hits remain. Never calculate the full mean by silently dropping censored observations. A standard error calculated from all hit times assumes independent trajectories; a fixed-horizon sample can use its restricted mean, but not claim a full MFPT when censored.
4. Compare quadrature, trajectory and grid. Discrepancies come from Monte Carlo sampling, missed crossings, finite step and finite grid. Verify their reduction before interpreting differences as physics. Grid endpoints have a finite-spacing boundary offset.
5. 1/MFPT is not automatically a phenomenological rate, and first hit of the barrier is not identical to reaching a product basin. Initial preparation matters.

## Facilitation

Ask students to predict each outcome before running it. Assign pairs to different barrier heights and D values; compare plots. Reserve five minutes to discuss molecular coordinates and hidden variables. An advanced challenge derives the variable-D drift and committor. For 90 minutes, use precomputed simulations or demonstrate them while students complete the derivation. Student tasks use TODO markers and include executable scaffolds so a fresh notebook runs; tasks explicitly ask students to edit/extend working code.

## Numerical and statistical notes

Euler–Maruyama uses Gaussian endpoint increments. Reflecting overshoots via a+|x−a| approximates a reflecting wall; an endpoint-only absorption misses bridge crossings. Check dt, dt/2 and dt/4 and report uncertainty. Never claim convergence just because curves look close. The teaching code is for moderate barriers and stable steps; extreme barriers may overflow exponentials and require log-domain quadrature. For a bounded domain with F=0, grid and continuum results agree in the limit h→0, not exactly for finite h.


## Foundations notebook: teaching notes and worked answers

Allow 60–90 additional minutes if this material is new. It can also be assigned as preparation. The foundations notebook has executable numerical experiments, prediction questions and a derivation exercise. The 150-minute schedule in the README describes the later session, not the additional foundations time.

1. For independent draws, sample SD approaches the population SD (one for a standard normal), while the standard error of the mean decreases approximately as $N^{-1/2}$. The particular sample mean fluctuates across repetitions; a standard error is not a deterministic error bound.
2. For free Brownian motion, the mean remains at the initial point and variance is $2Dt$. If the discrete increment instead has amplitude proportional to $\Delta t$, the variance after $t/\Delta t$ steps is proportional to $t\Delta t$ and tends to zero. Independent paths give independent endpoint samples; consecutive observations on one path are correlated.
3. With $\lambda=D\beta\kappa$, take $f(x)=x$: $\mathcal Lf=-\lambda x$, so $m'=-\lambda m$ and $m=x_0e^{-\lambda t}$. With $f=x^2$, $\mathcal Lf=-2\lambda x^2+2D$. Thus $M_2'=-2\lambda M_2+2D$, and for fixed $x_0$, $M_2=x_0^2e^{-2\lambda t}+(D/\lambda)(1-e^{-2\lambda t})$. Subtract $m^2$ to get the notebook's variance. The Euler scheme has finite-step bias; exact curves are continuous-time references, not exact predictions for the discretization.
4. In the exponential example, the observed-events-only average selects early arrivals. The restricted mean is $\int_0^T e^{-t/4}dt=4(1-e^{-T/4})$, not 4. The common horizon gives no information about the survival tail beyond it without additional data or a model.

The original root-level `01_student.ipynb` is an independently saved notebook; the course route in the README uses the notebooks under `notebooks/`. Avoid replacing a student's existing outputs or answers when distributing updates.

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

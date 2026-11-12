# Foundations: probability, Brownian motion and stochastic kinetics

This chapter assumes calculus and elementary Python, but no previous course in stochastic processes. Read it before `derivations.md`; the accompanying `notebooks/00_foundations.ipynb` lets you check its statements computationally. Symbols are defined where first used. The calculations below are derived under the stated smoothness and integrability assumptions; a measure-theoretic construction of Brownian motion is beyond the scope of this lesson.

## 1. Random variables, averages, and densities

A **random variable** $X$ associates a numerical result with each outcome of an experiment. Repeating the preparation produces a distribution of possible values. Its mean $\mathbb E[X]$ describes the ensemble average and $\operatorname{Var}(X)=\mathbb E[(X-\mathbb E[X])^2]$ its spread. The standard deviation is the square root of the variance and has the same units as $X$.

For a continuous random variable, a probability density $p(x)$ satisfies

$$\Pr(c<X<d)=\int_c^d p(x)\,dx,\qquad \int p(x)\,dx=1.$$

A density has units of inverse $x$; $p(x)$ is not the probability of observing exactly $x$. For independent samples $X_1,\dots,X_N$, the sample mean $\bar X=N^{-1}\sum_iX_i$ estimates $\mathbb E[X]$. Its standard error can be estimated as $s/\sqrt N$, where $s^2=\sum_i(X_i-\bar X)^2/(N-1)$. The **standard deviation** describes individual variation; the **standard error** describes uncertainty in the estimated mean. Adjacent frames of one trajectory are generally correlated and cannot automatically be treated as independent samples.

A conditional distribution, e.g. $p(x,t\mid x_0,0)$, describes position $x$ at time $t$ given the preparation $X_0=x_0$. Conditioning specifies the experiment; different initial preparations need not have the same kinetics.

## 2. What is a stochastic process?

A process $(X_t)_{t\ge0}$ is a family of random variables indexed by time. A single observed sequence $t\mapsto X_t(\omega)$ is a **sample path**; many independently prepared paths form an ensemble. The time-dependent density $p(x,t)$ describes that ensemble. Deterministic and stochastic models can both make precise predictions, but the latter predict distributions of outcomes.

A time-homogeneous **Markov process** has the property that, conditional on its present state, the future is independent of its earlier history. For transition densities this implies the Chapman–Kolmogorov relation

$$p(x,t+s\mid x_0,0)=\int p(x,t\mid y,0)p(y,s\mid x_0,0)\,dy.$$

Markovianity depends on the chosen state: the full position–momentum dynamics can be Markovian while position alone retains memory. In particular, an arbitrary molecular collective variable need not be Markovian just because the atomistic dynamics is.

## 3. Brownian motion and diffusion

Standard Brownian motion $W_t$ starts at zero, has continuous paths and independent increments, with $W_{t+\Delta t}-W_t\sim \mathcal N(0,\Delta t)$. This means a numerical increment is $\sqrt{\Delta t}\,\xi$ for $\xi\sim\mathcal N(0,1)$. Brownian paths are almost surely nowhere differentiable. Writing $\dot W$ as white noise is useful shorthand but it is not an ordinary time derivative.

The free-diffusion model

$$dX_t=\sqrt{2D}\,dW_t,\qquad X_0=x_0,$$

has solution $X_t=x_0+\sqrt{2D}W_t$ and therefore

$$\mathbb E[X_t]=x_0,\qquad \operatorname{Var}(X_t)=2Dt,\qquad [D]=\frac{[x]^2}{[t]}.$$

For $t>0$ its density is $p(x,t\mid x_0,0)=(4\pi Dt)^{-1/2}\exp[-(x-x_0)^2/(4Dt)]$. Substitution verifies the diffusion equation $\partial_t p=D\partial_x^2p$. The mean stays at the origin while individual realizations wander. In particular the typical displacement is proportional to $\sqrt{t}$, whereas its variance is proportional to $t$.

## 4. Langevin dynamics and the thermal environment

For a particle with position $x$, velocity $v$, mass $m$, friction coefficient $\zeta$, and potential energy $U$, a standard equilibrium model is

$$dx=v\,dt,\qquad m\,dv=[-U'(x)-\zeta v]dt+\sqrt{2\zeta k_BT}\,dW_t.$$

The stochastic forcing and friction satisfy fluctuation–dissipation at temperature $T$. Under suitable timescale separation, velocity relaxes on $m/\zeta$ much faster than the position changes appreciably. Eliminating the fast velocity gives the **overdamped** model for constant mobility $\mu=1/\zeta$,

$$dX_t=-\mu U'(X_t)dt+\sqrt{2\mu k_BT}\,dW_t
=-D\beta U'(X_t)dt+\sqrt{2D}\,dW_t,$$

where $D=\mu k_BT$ and $\beta=(k_BT)^{-1}$. This limit requires more than setting acceleration to zero in an arbitrary trajectory; it is a coarse-grained approximation. For a molecular coordinate, replacing $U$ by a free energy $F$ and using a Markov diffusion is an additional modeling assumption.

We choose length $\ell$, time $\ell^2/D_0$ and energy $k_BT$ as units. The resulting variables and equations in the notebook are dimensionless. The dimensionless double well $F(x)=H(x^2-1)^2$ has minima at $\pm1$ and a barrier of height $H$ at zero. Increasing $H$ in these units increases the barrier relative to thermal energy.

## 5. Itô increments and the generator

An Itô SDE is a definition through stochastic integrals, written informally as

$$dX_t=a(X_t)dt+\sqrt{2D(X_t)}\,dW_t.$$

In an Itô discretization the coefficients are evaluated at the beginning of each interval. For a smooth function $f$, Taylor expansion gives

$$\Delta f=f'\Delta X+\tfrac12 f''(\Delta X)^2+\cdots.$$

Since $\Delta W=O(\sqrt{\Delta t})$, the quadratic term contributes at order $\Delta t$: $\mathbb E[(\Delta W)^2]=\Delta t$. Terms of higher order vanish under the usual limiting assumptions. Itô's formula is

$$df=[a f'+D f'']dt+\sqrt{2D}f'\,dW.$$

Taking expectations (assuming the stochastic integral has mean zero) defines the **backward generator** $\mathcal L f=af'+Df''$ through

$$\frac{d}{dt}\mathbb E[f(X_t)]=\mathbb E[(\mathcal Lf)(X_t)].$$

This is the step that distinguishes stochastic calculus from simply applying the ordinary chain rule. The generator specifies evolution of observables, including quantities viewed as functions of the initial position.

## 6. From probability conservation to Fokker–Planck

Write $\mathbb E[f(X_t)]=\int f(x)p(x,t)dx$. For smooth $p,f$, and boundary conditions that eliminate the boundary terms,

$$\int (af'+Df'')p\,dx=\int f[-\partial_x(ap)+\partial_x^2(Dp)]\,dx.$$

Thus the **forward** equation is

$$\boxed{\partial_t p=-\partial_x(ap)+\partial_x^2(Dp)=-\partial_xJ},\qquad J=ap-\partial_x(Dp).$$

This is a continuity equation: the probability inside an interval changes through current across its boundaries. For constant $D$ and $a=-D\beta F'$,

$$J=-D(\beta F'p+p'),\qquad p_{\rm eq}=Z^{-1}e^{-\beta F}$$

has zero current when normalizable. A reflecting boundary has $J=0$; a regular absorbing endpoint has $p=0$. At equilibrium, $D$ changes the time scale but not this normalized distribution. Stationarity need not mean zero current in every nonequilibrium or periodic problem; here we focus on equilibrium detailed balance.

**Variable $D$ and conventions.** For an equilibrium Itô diffusion with specified $D(x)$ and $p_{\rm eq}\propto e^{-\beta F}$, the drift is $a=D'-D\beta F'$. This follows directly by setting $J=0$. Omitting $D'$ changes the model. Different stochastic integral conventions require an appropriately transformed drift.

## 7. First passage, survival and mean time

For a process started at $x$, define $\tau_b=\inf\{t\ge0:X_t\ge b\}$. Its survival function is $S_x(t)=\Pr_x(\tau_b>t)$. If $\mathbb E_x\tau_b<\infty$ then

$$\mathbb E_x\tau_b=\int_0^\infty S_x(t)dt.$$

To see this, write $\tau_b=\int_0^\infty\mathbf1_{\{\tau_b>t\}}dt$ and exchange expectation and integration. If trajectories are stopped at a finite horizon $T$, an unobserved event is **right-censored**. The average $\langle\min(\tau_b,T)\rangle$ estimates the *restricted* mean $\int_0^T S(t)dt$, not the full MFPT; averaging only observed events is generally biased. A fixed administrative censoring horizon allows an empirical survival curve strictly before the horizon. At the horizon, a censored record does not tell whether $\tau=T$ or $\tau>T$.

With a reflecting left boundary $a$ and absorbing target $b$, first-step analysis gives $\mathcal L T=-1$, with $T'(a)=0$ and $T(b)=0$. `derivations.md` works through the integration. Hitting $b$ is an event definition; hitting a barrier top and reaching a product basin are different. $1/\mathbb E\tau$ is not automatically a phenomenological rate constant without suitable preparation, metastability and recrossing assumptions.

## 8. Reading list

Use `derivations.md` next for the integrating factor, committor, variable diffusion and grid generator. `references.md` gives sources and attribution. The notebook contains explicit numerical experiments and uncertainty checks. For formal stochastic calculus see Gardiner; for molecular reaction coordinates and position-dependent diffusion see Hummer (2005), both in the reference list.

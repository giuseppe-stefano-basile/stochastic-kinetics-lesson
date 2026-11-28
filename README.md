# Stochastic kinetics: from Langevin trajectories to Fokker–Planck

A hands-on, mathematically grounded introduction for advanced undergraduates and beginning graduate students. Starting from a single overdamped particle, you will derive the forward and backward equations, calculate first-passage times, and test predictions against simulation and a discrete generator.

## Learning goals
- Connect trajectories, transition probabilities and probability density.
- Derive the Itô generator, Fokker–Planck equation and probability current.
- Distinguish equilibrium from kinetics and explain boundary conditions.
- Compute and validate mean first-passage times (MFPTs) and committors.
- Recognize when a molecular collective variable may not be Markovian.

## Quick start

Use Python 3.10 or newer. From the repository root on Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[teaching]'
jupyter lab
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1` after creating the environment with `py -m venv .venv`. In Jupyter, open the notebooks under `notebooks/` and use the Python kernel from this environment. Once installed, the examples run locally without external datasets. Run `python -m pytest` from the repository root to check the numerical code.

## Lesson path

| Material | What to do |
| --- | --- |
| [Foundations notes](notes/00_foundations.ipynb) and [foundations lab](notebooks/00_foundations.ipynb) | Build intuition for distributions, Brownian paths and ensemble variance. |
| [Student notebook](notebooks/01_student.ipynb) | Model a double well, estimate survival and compare first-passage methods. |
| [Kinetics investigation](notebooks/03_kinetics_investigation.ipynb) | Match the model and boundaries across methods, then change the diffusion, barrier or target. |

The [derivations](notes/derivations.ipynb), [instructor guide](notes/instructor.ipynb), and [references](notes/references.ipynb) are Jupyter notebooks with sectioned Markdown cells and rendered equations. The [worked notebook](notebooks/02_instructor.ipynb) contains executable solutions. The [printable instructor notes](notes/Stochastic_Kinetics_Instructor_Notes.pdf) are also available as [LaTeX source](notes/Stochastic_Kinetics_Instructor_Notes.tex).

The model uses dimensionless units with default `beta=1` and `D=1`. The trajectory code reflects overshoots and checks absorption at sampled endpoints. The grid generator approximates the continuum with a finite-state jump process. Compare step sizes and grids, and account separately for sampling uncertainty and any censored observations.

## Attribution and license

Inspired by Luca Sagresti's [HandsOnSession_CMLS](https://github.com/lianctrl/HandsOnSession_CMLS), particularly its Langevin, MSM and Fokker–Planck tutorial. This project implements its own model and code and does not redistribute files from that repository. Consult the original repository's license before reusing its material. Original material here uses the [MIT license](LICENSE). See the [figure credits](slides/SNS_BUILD.md#figure-credits) for reproduced figures and their license.


## Beamer slides

Open the [compiled Beamer slides](slides/stochastic_kinetics.pdf) or edit [the LaTeX source](slides/stochastic_kinetics.tex). The main deck opens with molecular examples and gives explicit notebook handoffs. Additional derivations and extensions follow the **Backup slides** divider.

Build with XeLaTeX and a TeX distribution containing Beamer, TikZ, PGFPlots and fontspec:

```bash
cd slides
xelatex -halt-on-error stochastic_kinetics.tex
xelatex -halt-on-error stochastic_kinetics.tex
```

See [SNS build and figure credits](slides/SNS_BUILD.md) for fonts and attribution.

For the Beamer density comparison, `python slides/generate_density_data.py` regenerates the CSV data in the working directory. The slide source embeds the displayed coordinates, so XeLaTeX can build the deck without running Python.

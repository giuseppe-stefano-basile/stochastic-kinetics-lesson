# Stochastic kinetics: from Langevin trajectories to Fokker–Planck

A hands-on, mathematically grounded introduction for advanced undergraduates and beginning graduate students. Starting from a single overdamped particle, you will derive the forward and backward equations, calculate first-passage times, and test predictions against simulation and a discrete generator.

## Learning goals
- Connect trajectories, transition probabilities and probability density.
- Derive the Itô generator, Fokker–Planck equation and probability current.
- Distinguish equilibrium from kinetics and explain boundary conditions.
- Compute and validate mean first-passage times (MFPTs) and committors.
- Recognize when a molecular collective variable may not be Markovian.

## Quick start

Python 3.10+ recommended. In a terminal:

```bash
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
python -m pip install -e '.[teaching]'
jupyter lab
```

Start with `notes/00_foundations.md` and `notebooks/00_foundations.ipynb`, then open `notebooks/01_student.ipynb`. Alternatively upload the entire repository to a notebook environment and install with `pip install -e .`; no network or external dataset is needed after installation. The student notebook contains tasks and runnable starter cells; `02_instructor.ipynb` gives worked code and results. Read `notes/derivations.md` alongside them. Execute cells in order. Plotting requires a notebook kernel. Run `python -m pytest` for automated checks.



Students new to probability should complete the foundations chapter and notebook before class, if possible.

## Contents

- `notes/00_foundations.md`: a beginner introduction to probability, Brownian motion and the Fokker–Planck equation.
- `notebooks/00_foundations.ipynb`: worked numerical introduction and guided exercises.
- `notebooks/01_student.ipynb`: scaffold, predictions and exercises.
- `notebooks/02_instructor.ipynb`: executable solutions and diagnostics.
- `notes/derivations.md`: detailed mathematical development and assumptions.
- `notes/instructor.md`: timings, expected answers, choices and pitfalls.
- `notes/references.md`: references and repository attribution.
- `src/stochkin/`: short transparent numerical implementations.
- `tests/`: analytic and numerical validation.

The model uses dimensionless units: beta=1, length scale=1, D=1 by default. A simulation step observes the process only at discrete times; reflection by folding and absorption by endpoint crossing approximate continuous boundaries. Check convergence in the step size. The grid generator is a continuous-time Markov-chain approximation and its boundary is a grid node, not an exact continuum boundary at finite spacing.

## Attribution and license

Inspired by Luca Sagresti's [HandsOnSession_CMLS](https://github.com/lianctrl/HandsOnSession_CMLS), particularly its Langevin, MSM and Fokker–Planck tutorial. This project implements its own model and code and does not redistribute files from that repository. Consult the original repository's license before reusing its material. Original material here is offered under the MIT license in `LICENSE`.


## Beamer slides

The editable presentation is `slides/stochastic_kinetics.tex`. Compile it with a standard TeX Live installation containing Beamer, TikZ and PGFPlots:

```bash
cd slides
pdflatex stochastic_kinetics.tex
pdflatex stochastic_kinetics.tex
```

The slides include the introductory material and the main hands-on lesson. Plan an additional 60–90 minutes for students who are new to stochastic processes, or assign the foundations as preparation.

### Guided kinetics investigation

Open `notebooks/03_kinetics_investigation.ipynb` for a single double-well experiment connecting the backward equation, simulated first passages, grid refinement, censoring and parameter changes. For the Beamer density comparison, `python slides/generate_density_data.py` regenerates the CSV data in the current directory (NumPy and SciPy required). The slide source embeds the displayed coordinates and builds without Python. The enlarged slide deck is a full teaching resource; use the introduction as preparation and select advanced slides according to the available class time.

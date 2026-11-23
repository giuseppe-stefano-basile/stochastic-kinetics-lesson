# Build the Beamer deck

From this directory, run XeLaTeX twice:

```bash
xelatex -halt-on-error stochastic_kinetics.tex
xelatex -halt-on-error stochastic_kinetics.tex
```

Install the SNS title and body fonts (Rotis Semi Serif and Roboto Condensed) to use the institutional pair. The theme uses Arial when the pair is unavailable and Liberation Sans when Arial is also unavailable. The PDF supplied here uses Liberation Sans.

The Beamer theme follows the public [SNS visual identity guide](https://www.sns.it/it/guida/identita-visiva). The official PowerPoint/Keynote templates are distributed separately. If SNS Communications has approved logo use for this presentation, place the approved horizontal PDF next to the TeX source as `sns-logo-approved.pdf` and rebuild.

## Figure credits

Figures 1 and 6 come from Sagresti, Peri, Ceccarelli and Brancato, *J. Chem. Theory Comput.* **18** (2022), 3164–3173, [DOI 10.1021/acs.jctc.2c00181](https://doi.org/10.1021/acs.jctc.2c00181). The image panels were extracted unchanged from the published PDF. They are © 2022 the authors, published by ACS, under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

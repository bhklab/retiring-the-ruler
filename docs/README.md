# Retiring the Ruler

**Authors:** [Katy Scott](https://github.com/strixy16), Caryn Geady, Kaitlyn Kobayashi

**Contact:** [bhklab.katyscott@gmail.com](mailto:bhklab.katyscott@gmail.com)

**Description:** Simulation code accompanying the Retiring the Ruler: Toward a Volumetric, Multi-Lesion, Longitudinal Lens on Tumor Response paper.

--------------------------------------

[![pixi-badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json&style=flat-square)](https://github.com/prefix-dev/pixi)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json&style=flat-square)](https://github.com/astral-sh/ruff)
[![Built with Material for MkDocs](https://img.shields.io/badge/mkdocs--material-gray?logo=materialformkdocs&style=flat-square)](https://github.com/squidfunk/mkdocs-material)

![GitHub last commit](https://img.shields.io/github/last-commit/bhklab/retiring-the-ruler?style=flat-square)
![GitHub issues](https://img.shields.io/github/issues/bhklab/retiring-the-ruler?style=flat-square)
![GitHub pull requests](https://img.shields.io/github/issues-pr/bhklab/retiring-the-ruler?style=flat-square)
![GitHub contributors](https://img.shields.io/github/contributors/bhklab/retiring-the-ruler?style=flat-square)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/bhklab/retiring-the-ruler?style=flat-square)

## Set Up

### Prerequisites

Pixi is required to run this project.
If you haven't installed it yet, [follow these instructions](https://pixi.sh/latest/)

### Installation

1. Clone this repository to your local machine
2. Navigate to the project directory
3. Set up the environment using Pixi:

```bash
pixi install
```

## Algorithm
### SYNTHETIC LESION GENERATION
1. Artificially generate a dataset of N patients (input arg: N).
2. For each patient, generate a random number of lesions (range: 1-30, poisson distribution with µ = 5).
3. For each lesion, generate a random diameter (select from a list of possible diameters, with replacement (using original_shape_Maximum2DDiameterSlice)).
4. For each lesion, generate a diameter change (select from a Gaussian distribution, with fractional changes from -1 to 1).
5. For each lesion, assign a location tag (e.g. liver, lung, bone, brain, again from a list of possible locations, with replacement).

### LESION SELECTION BIAS ASSESSMENT

6. For each patient, select a random number of lesions to be target lesions (up to 5 per patient, with a maximum of 2 per assigned location).
7. From the selected target lesions, assess the RECIST response (SLD = 100% decrease = CR, SLD > 30% decrease = PR, SLD > 20% increase = PD, otherwise SD).
8. For each patient, assess the RECIST response across all lesions (as above).
9. Compare the RECIST response across all lesions to the RECIST response across target lesions only.

### DIAMETER VERSUS VOLUME ASSESSMENT

10. For each diameter, calculate the corresponding volume (assuming spherical shape; this will represent the lower end of volumetric measurement).
11. Inject inter-observer variability by adding a random error to the volume calculation (select from a Gaussian distribution, with fractional diameter changes from -0.113 to 0.113; REF: https://pmc.ncbi.nlm.nih.gov/articles/PMC3423763/).
12. For each lesion, calculate the volume using another diameter measurement (i.e., original_shape_Maximum3DDiameter).


## Documentation

Click [here](https://bhklab.github.io/retiring-the-ruler) to view the full documentation.

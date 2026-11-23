# Fitness-cost figures

This code reproduces the quantitative figures in the manuscript and Supplementary Information. The experimental data are included in `experiment_data.xlsx`.

## Run

Python 3 is required. From this directory, install the dependencies and run:

```bash
python -m pip install numpy matplotlib pandas scipy pillow openpyxl
python build_all.py
```

The generated PNG files will be saved in `output/`.

To regenerate only one group of figures, run its module directly. For example:

```bash
python -m plotting.phir_gex
```

## Files

- `build_all.py`: generates all figures.
- `model.py`: model parameters and equations.
- `experiment_data.xlsx`: experimental input data.
- `plotting/`: data analysis, shared helpers, and figure modules.

## Notes

The data were digitized from published figures and may contain small extraction errors. These uncertainties do not affect the main conclusions.

For PGK1 and PDC1 in Fig. S4C of Kafri et al. (2016), overlapping data points could not be reliably extracted from the figure, so the published fit-line slopes were used.

# HERMES: Health, Epidemic and Economic R-based Microsimulation Engine for Individual-based Simulations

HERMES is a lightweight **individual-based modelling (IBM) framework** for the
simulation of infectious disease transmission in structured populations.
The framework explicitly represents **heterogeneous individuals** embedded in
households, schools, workplaces, and the wider community, and is designed to
support **research, teaching, and reproducible computational modelling**.

HERMES is intentionally developed as a collection of **plain R scripts** that are
transparent, modular, and easily adaptable, lowering the barrier for method
development, inspection, and extension.

---
## Installation (using R)

Download all R files using `download.file()`. For more info, please consult
https://github.com/lwillem/ibm_hermes/blob/main/install.R



## Installation (using GitHub)

Clone or download the repository `https://github.com/lwillem/ibm_hermes.git`.


---

## Quick start

Set your working directory to the project root and run the main workbench script:

```r
source("main.R")
```

This will:
- load all required functions and dependencies,
- initialise the default model parameters,
- execute an example simulation,
- optionally perform an internal regression test.

Alternatively, the model kernel can be invoked directly:

```r
source("lib/ibm_core.R")

params <- get_default_parameters()
out <- run_ibm(params)
```


---

## To use with modified parameters

Model behaviour is fully controlled via the `params` object.  
All available parameters can be inspected using:

```r
print_model_parameters()
```

For example, to change the number of initially infected individuals and enable
comparison with a baseline scenario:

```r
# load all helper functions
source("lib/ibm_core.R")

# retrieve default parameters
params <- get_default_parameters()

# modify selected parameters
params$num_infected_seeds <- 10     # increase number of initial infections
params$bool_add_baseline  <- TRUE   # enable baseline comparison

# run the model
out <- run_ibm(params)
```

---

## Model structure

The framework is organised into clearly separated components, each with a
well-defined responsibility:

| File | Responsibility |
|----|----|
| `main.R` | Workbench script for running model scenarios |
| `lib/ibm_core.R` | Core simulation kernel and internal helper functions |
| `lib/ibm_parameters.R` | Definition, inspection, and comparison of model parameters |
| `lib/ibm_population.R` | Synthetic population generation |
| `lib/ibm_plot.R` | Visualisation and plotting utilities |
| `lib/ibm_test.R` | Regression testing and reference management |

This modular structure facilitates code readability, reuse, and systematic
model extension.

---

## Model output

Simulation output is written to the `output` directory, in a subdirectory whose
name is specified via the `params` object. By default, output is stored in `output/ibm_flu`

Results are saved as RDS files and include:

- **`health_states.rds`**  
  A list containing:
  - `log_health`: a data frame with population-level health state proportions
    over time (columns represent health states; rows represent time steps),
  - `params`: the full set of model parameters used for the simulation.

- **`pop_matrix.rds`**  
  A data frame containing individual-level population data at the end of the
  simulation, with one row per individual and individual attributes stored in
  columns.

---

## Extending the model

To extend HERMES in a consistent and reproducible manner:

- Add or modify parameters in `get_default_parameters()`
- Extend individual attributes in `create_population_matrix()`
- Extend model logic within `run_ibm()`
- Introduce new internal helper functions (e.g. transmission or mortality logic)
  in `ibm_core.R`
- Maintain `run_ibm()` as the **single model entry point**

This design ensures that regression testing remains informative and that
structural changes are explicit.

---

## Reproducibility

- All stochastic components are controlled via `params$rng_seed`
- A lightweight internal regression test is available:

```r
run_ibm_regression_test()
```

- To intentionally update the stored reference output:

```r
update_ibm_reference()
```

---

## Author

**Lander Willem**  
University of Antwerp, Belgium

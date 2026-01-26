# HERMES: Health, Epidemic and Economic R-based Microsimulation Engine for individual-based Simulations

HERMES is a lightweight **individual-based modelling (IBM) framework** for
simulating infectious disease transmission in structured populations.
It supports households, schools, workplaces, and community transmission,
and is designed for **research, teaching, and reproducible modelling**.

The framework is intentionally designed to work both:
- as a collection of **plain R scripts**, and
- as a **minimal R package**, without forcing package workflows.

---

## Install

Clone or download the repository. Installation or compilation is not needed at this stage.

## Quick start

Set your working directory to the project root, you can run the main file with:

```r
source("main.R")
```

or:

```r
source('lib/ibm_core.R')

params <- get_default_parameters()
out <- run_ibm(params)
```

This will:
- load all required functions,
- initialise default parameters,
- run an example simulation,
- optionally perform a regression test.

---

## Quick start (with changed parameters)

If you want to change model parameters, you can do so by changing the content of the `params` variable. You can explore the options by printing all options with `print_model_parameters()`. For example, when changing the number of initial infections to 10:

```r
# load all helper functions
source('lib/ibm_core.R')

# get default parameters
params <- get_default_parameters()

# change some parameters
params$num_infected_seeds <- 10  # change infected cases to 10
params$bool_add_baseline  <- TRUE # option to show the difference with baseline

# run model
out <- run_ibm(params)
```

---

## Model structure

The framework is organised into clearly separated components:

| File | Responsibility |
|----|----|
| `main.R` | Workbench script for running scenarios |
| `ibm_core.R` | Core model kernel and internal helpers |
| `ibm_parameters.R` | Parameter definitions and exploration |
| `ibm_population.R` | Synthetic population generation |
| `ibm_plot.R` | Plotting utilities |
| `ibm_test.R` | Regression testing |

This separation makes the code easy to read, reuse, and extend.

---

## Model output

Model output is stored in the `output` directory in a subdirectory with a name included in the `params` variable. The current default output directory is `output/ibm_flu`. The output is saved as rds files and includes:

- health_states.rds A list containing a data.frame `log_health` with the health states in the columns and the time steps in the rows and list called `params` containing all parameters used to obtain this output.
- pop_matrix.rds The population details at the end of the simulation, with one row per person and the person details in the columns. 


## Extending the model

To adopt HERMES safely:

- Add or modify parameters in `get_default_parameters()`
- Extend individual attributes in `create_population_matrix()`
- Extend the model logic in `run_ibm()`
- Add new internal helpers (e.g. transmission or mortality logic) in
  `ibm_core.R`
- Keep `run_ibm()` as the **single model entry point**

This design ensures regression tests remain meaningful.

---

## Reproducibility

- All stochastic behaviour is controlled via `params$rng_seed`
- A lightweight regression test is provided:

```r
run_ibm_regression_test()
```

- To update the reference output intentionally:

```r
update_ibm_reference()
```

---

## Author

**Lander Willem**  
University of Antwerp

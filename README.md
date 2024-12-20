# Simulation Codes: Unit Averaging for Heterogeneous Panels

## Overview


This repository contains the simulation codes used in the paper "Unit Averaging for Heterogeneous Panels" by Christian Brownlees and Vladislav Morozov. The codes are written in MATLAB and can be run to reproduce the simulation results presented in the [paper](https://arxiv.org/abs/2210.14205) and the [Online Appendix](https://vladislav-morozov.github.io/assets/files/1_unitAveragingSupplement.pdf).

## Contents


 
  - [About: What Is Unit Averaging](#about-what-is-unit-averaging)
  - [Results: Unit Averaging Improves Efficiency](#results-unit-averaging-improves-efficiency)
  - [Software Requirements](#software-requirements)
  - [How to Replicate the Simulations](#how-to-replicate-the-simulations)
  - [Folder Structure](#folder-structure)
  - [License](#license)
  - [External Functions Used](#external-functions-used)
  - [Contact](#contact)



## About: What Is Unit Averaging

### What Is Unit Averaging?

In short, **unit averaging** is a simple approach to efficiently estimate unit-specific parameters in heterogeneous panel data settings. It can easily and automatically accommodate both linear and nonlinear models.



### When to Use Unit Averaging?
 
Use unit averaging when the following two conditions hold:
1. You need to estimate a parameter of a specific business or economic unit, such as countries, customers, or firms.
2. You have a panel of similar units at your disposal. 
   
Unit averaging can be used both in causal and forecasting frameworks. Any unit-specific quantity may be treated as the parameter of interest.

In any setting of interest, the panel is likely to be heterogeneous in the the sense that every unit may have a different value of the parameter. Unit averaging should be used in such settings to efficiently use the information in this panel while taking heterogeneity into account.

### Example Settings

  - Forecasting GDP for France given a panel of European countries.
  - Estimating the production function of the plastic producing sector in the US given a panel of chemical producers.
  - Estimating the change in the probability that a given customer buys a product after a promotion, given a panel of other customers.



### How Does Unit Averaging Work?

Unit averaging proceeds as follows:

1. Estimate the parameter of interest separately for each unit in the panel.
2. Combine the unit-specific parameter estimates using a linear combination.

We propose an optimal weight scheme that minimizes a mean squared error (MSE) criterion.

 ## Results: Unit Averaging Improves Efficiency

Unit averaging with our optimal weights leads to more precise estimator compared both to no averaging and other weight specifications.  It does so by optimally striking the bias-variance trade-off.

For example, the image below shows that our optimal unit averaging reduce the MSE in the whole parameter space, comparing to no averaging. It does so by optimally trading-off bias and variance. Optimal weights do not suffer from breakdowns in performance, unlike equal weights (the mean group estimator).

![Simplified simulation results: out unit averaging approach leads to lower mean squared error](results/figures/animated_simplified_results.gif?raw=true)
 

## Software Requirements
- MATLAB (tested with MATLAB 2022b and 2024b). 
- Required toolboxes: Parallel Computation Toolbox.
- Required File Exchange files: `table2latex` and `tight_subplot`. Supplied with the replication code. 


## How to Replicate the Simulations
1. Clone the repository to your local machine.
2. Open MATLAB and set the current folder to the repository directory.
3. Run the `main.m` file to execute the full simulations.

## Folder Structure

- `main.m`: The main script to run the simulations.
- `results/`:
  - `figures/`: Contains the exported plots.
  - `simulations/`: Contains the `.mat` files with simulation results.
  - `tables/`: Contains the exported tables.
- `src/`: Contains source code files organized into subdirectories based on their goals.
  - `averaging_functions/`: Unit averaging weight scheme functions.
  - `configs/`: Config scripts with simulation and plotting parameters.
  - `utilities/`: Assorted functions:
    - `averaging/`: auxiliary functions for unit averaging.
    - `results_plots_tables`: functions for processing and exporting results.
    - `sampling`: functions for drawing and resampling data.
- `scripts`: Contains the main simulation and exporting scripts
  - `plotting/`: Scripts for exporting plots.
  - `simulations/`: Simulation scripts
  - `tables/`: Scripts for exporting tables.

> [!NOTE]
> The replication folder does not contain the rather large final simulation files, figures, and tables. These files can be generated by running `main.m`


## License
This code is provided under the MIT License. You can use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the software. You must include the original copyright notice and this permission notice in all copies or substantial portions of the software. The software is provided "as is", without any warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement.

## External Functions Used

We have grateful to the creators of the following functions on the Matlab File Exchange:
1. `latex2table` by Víctor Martínez-Cagigal (https://www.mathworks.com/matlabcentral/fileexchange/69063-matlab-table-to-latex-conversor);
2. `tight_subplot` by Pekka Kumpulainen (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w).

## Contact

For questions or feedback, please contact the authors:
- [Christian Brownlees](https://github.com/ctbrownlees/)
- [Vladislav Morozov](https://github.com/vladislav-morozov)
 
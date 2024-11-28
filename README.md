# Simulation Codes for Unit Averaging for Heterogeneous Panels

## Overview


This repository contains the simulation codes used in the paper "Unit Averaging for Heterogeneous Panels" by Christian Brownlees and Vladislav Morozov. The codes are written in MATLAB and can be run to reproduce the simulation results presented in the [paper](https://arxiv.org/abs/2210.14205) and the [Online Appendix](https://vladislav-morozov.github.io/assets/files/1_unitAveragingSupplement.pdf).

## Contents



 
- [What is Unit Averaging](#what-is-unit-averaging)
- [Software Requirements](#software-requirements)
- [How to Run the Simulations](#how-to-run-the-simulations)
- [Directory Structure](#directory-structure)
- [License](#license)
- [Contact](#contact)


## About Unit Averaging

In short, **unit averaging** is a simple approach to efficiently estimation unit-specific parameters using heterogeneous panel data models. It can easily accommodate both linear and nonlinear models.

Parameters of interest include

In essence, unit averaging proceeds as follows:
1. Individual
2. Averaging

The weights are determined by minimizing a mean squared error (MSE) criterion we derive. 

 ## Results: Unit Averaging Improves Efficiency

Optimally trades off bias and variance

Example image

## Software Requirements
- MATLAB (tested with MATLAB 2022b and 2024b). 
- Required toolboxes: Parallel Computation Toolbox.


## How to Replicate the Simulations
1. Clone the repository to your local machine.
2. Open MATLAB and set the current folder to the repository directory.
3. Run the `main.m` file to execute the simulations.

## Directory Structure

- `main.m`: The main script to run the simulations.
- `results/`: 
- `src/`: Contains the source code files. 
  - ` /`: Subdirectory for  
- `scripts`: 

## License
This code is provided under the MIT License. 

## Contact

For questions or feedback, please contact the authors:
- [Christian Brownlees](https://github.com/ctbrownlees/)
- [Vladislav Morozov](https://github.com/vladislav-morozov)
 
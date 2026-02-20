![Made with MATLAB](https://img.shields.io/badge/Made%20with-MATLAB-orange)
![Research Project](https://img.shields.io/badge/Project-Research-blue)
![License MIT](https://img.shields.io/badge/License-MIT-green)

# Changes in EMG–angle relationships following reverse total shoulder arthroplasty

This repository contains the MATLAB code used for the analyses of Changes in EMG–angle relationships following reverse total shoulder arthroplasty. 

## Overview

This project investigates changes in the relationship between thoracohumeral (TH) elevation and shoulder muscle activation before and after reverse total shoulder arthroplasty (rTSA), compared with asymptomatic subjects.

Surface EMG and kinematic data recorded during standardised functional tasks were used to quantify EMG–angle coupling through linear regression and variability analysis (RMSD).

## Repository content

- EMG and kinematic preprocessing  
- EMG–angle regression analysis  
- Variability (RMSD) computation  
- Statistical analyses and figure generation  

## Data structure (not publicly available)

The code relies on preprocessed `.mat` files that are **not included** in this repository due to ethical constraints.

The main data files used in the analyses are:

### EMG data
These files contain **EMG signals** for the combined functional tasks.

- `combined_functional_data_ASYMPTO.mat`
- `combined_functional_data_PREOP.mat`
- `combined_functional_data_POSTOP.mat`

### Kinematic-only data
These files contain **thoracohumeral elevation angle** for the combined functional tasks.

- `all_angles_statistics_Asymptomatic.mat`
- `all_angles_statistics_Pre_operatoire.mat`
- `all_angles_statistics_Post_operatoire.mat`

Additional files:
- `Cross_analysis_EMGxHT.m`: Main script (EMG–angle coupling analysis)
- `plotCombinedEMGPerSubjectWithSPM1D.m`: EMG analysis and mean cycle for the combined functional tasks
- `cycle_selections_combined_functional_data.mat`: selected movement cycles used for averaging the mean cycle

## Methods

- EMG: band-pass filtering (15–475 Hz), rectification, RMS smoothing (250 ms), submaximal-task normalisation, time-normalised to 100% of the movement cycle  
- Kinematics: TH elevation computed using ISB-recommended Euler angles  
- Analysis: linear regression slopes and RMSD of EMG–angle relationships  
- Statistics: ANOVA and post-hoc tests

## Requirements

- MATLAB R2024b  
- Signal Processing and Statistics Toolboxes

## Data availability

Data can be made available upon reasonable request to the corresponding author.

## Contact

**Florent Moissenet**  
florent.moissenet@unige.ch

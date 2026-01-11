# Changes in EMG–angle relationships following reverse total shoulder arthroplasty

This repository contains the MATLAB code used for the analyses presented in the abstract:

**Francalanci H., Holzer N., Cherni Y., Moissenet F.**  
*Changes in EMG–angle relationships following reverse total shoulder arthroplasty*  
ISG Lyon 2026.

## Overview

This project investigates changes in the relationship between thoracohumeral (TH) elevation and shoulder muscle activation before and after reverse total shoulder arthroplasty (rTSA), compared with asymptomatic subjects.

Surface EMG and kinematic data recorded during standardised functional tasks were used to quantify EMG–angle coupling through linear regression and variability analysis (RMSD).

## Repository content

- EMG and kinematic preprocessing  
- EMG–angle regression analysis  
- Variability (RMSD) computation  
- Statistical analyses and figure generation  

## Methods (summary)

- EMG: band-pass filtering (15–475 Hz), rectification, RMS smoothing (250 ms), submaximal-task normalisation, time-normalised to 100% of the movement cycle  
- Kinematics: TH elevation computed using ISB-recommended Euler angles  
- Analysis: linear regression slopes and RMSD of EMG–angle relationships  
- Statistics: ANOVA and post-hoc tests (α = 0.05)

## Requirements

- MATLAB R2024b  
- Signal Processing and Statistics Toolboxes

## Data availability

Raw data are not publicly available due to ethical constraints. Data can be provided upon reasonable request to the corresponding author.

## Contact

**Florent Moissenet**  
florent.moissenet@unige.ch

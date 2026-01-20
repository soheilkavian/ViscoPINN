# ViscoPINN
MATLAB code for physics-informed viscosity prediction of multicomponent aqueous solutions

# ViscoPINN: Physics-Informed Viscosity Prediction Software

## Overview
**ViscoPINN** is a MATLAB-based software package for predicting temperature-dependent viscosity of multicomponent aqueous mixtures using a **Physics-Informed Neural Network (PINN)** framework. The software supports mixtures of up to **17 pure components** and enables both direct PINN inference and optional **Vogel–Fulcher–Tammann (VFT)** extrapolation toward the glass transition regime.

The package is distributed in two complementary forms:
1. A **standalone Windows executable (.exe)** for plug-and-play use without MATLAB  
2. A **MATLAB function-based interface** for scripting, automation, and reproducibility

---

## Package Contents

### Core Files
- **`ViscosityPredictor.m`**  
  Main MATLAB function implementing PINN-based viscosity prediction, including input normalization, network inference, optional VFT fitting, and output formatting.

- **`PINN.mat`**  
  Trained Physics-Informed Neural Network model required for all predictions.

- **`mixture_data.mat`**  
  Supporting mixture property data used internally by the model.

- **`Tg_fitted_VFT.mat`**  
  Pre-fitted parameters used for VFT-based extrapolation.

---

### Example and Output Files
- **`Example.m`**  
  Example MATLAB script demonstrating how to call `ViscosityPredictor`, define mixture composition, specify temperature range, and visualize results.

- **`MyResults.xlsx`**  
  Example output file generated when exporting predicted viscosity–temperature data.

---

### Standalone Application
- **`MyAppInstaller_web.exe`**  
  Standalone Windows executable providing a graphical user interface (GUI) for viscosity prediction.  
  This version does **not require MATLAB** and is intended for user-friendly, interactive use.

---

## Standalone GUI Application (No MATLAB Required)

The ViscoPINN GUI enables interactive viscosity prediction through a plug-and-play Windows interface.

### Interface Structure

#### Input Panel
- Define mixture composition using **mole fractions or mass fractions**
- Supports up to **17 pure components**
- Input fractions are internally normalized

#### Control Panel
- Specify temperature range:
  - Lower bound (°C)
  - Upper bound (°C)
  - Step size (°C)
- **VFT Fit option**:
  - **Unchecked (default):** Direct PINN inference  
    - Highest accuracy within the training range of **−20 °C to 35 °C**
  - **Checked:** PINN + VFT fitting  
    - Enables theoretical extrapolation toward the glass transition temperature
- Buttons:
  - **Run Prediction**
  - **Export to Excel**

#### Output Panel
- Displays predicted viscosity as a continuous function of temperature
- Output units: **centipoise (cP)**
- Enables rapid assessment of non-linear temperature and composition effects

---

## MATLAB Function Usage

Advanced users may call the prediction engine directly from MATLAB for scripting, batch processing, or integration with other models.

### Function Call
```matlab
finalResults = ViscosityPredictor( ...
    moleFractions, ...
    LowT, StepT, HighT, ...
    useVFT, ...
    outputFile );

## Citation

If you use this software, please cite:

Soheil Kavian, Arian Zarriz, Matthew J. Powell-Palm1, [Journal Name], 2026.
A Large-Scale Dataset and Physics-Informed Neural Network for Viscosity Prediction in Many-Component Aqueous and Organic Solutions.

(Full citation provided in the manuscript.)

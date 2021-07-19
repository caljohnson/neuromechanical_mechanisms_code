# MATLAB Code for two publications:

Carter L. Johnson (2020) "Neuromechanical Mechanisms of Locomotion in C. elegans" Ph.D. Dissertation, UC Davis. DOI: 10.13140/RG.2.2.24387.02082

Carter L. Johnson, Timothy J. Lewis, and Robert D. Guy (2021), "Neuromechanical Mechanisms of Gait Adaptation in C. elegans: Relative Roles of Neural and Mechanical Coupling", SIAM Journal on Applied Dynamical Systems. DOI: 10.1137/20M1346122

## How to Generate Each Figure

### Figures from 2021 Paper

### Figures from 2020 Dissertation


## Folders

* ./ - Main codes to run simulations, parameter sweeps, generate data and figures
* ./src - Source codes and functions that main codes run off of

## src codes

### Model Codes
* full_state_model_odes.m - function that generates ODEs for the forward locomotion system (n neuromechanical modules)
* coupled_oscillator_phase_difference_odes.m - function that generates phase-equation ODEs for the forward locomotion system (n oscillators)
* coup_currents_and_oscillator_coupling_fns.m / oscillator_coupling_fns.m - computes the coupling functions for the full coupled system from the single oscillator limit cycle and PRC
* single_oscillator_LC.m - computes the limit cycle oscillations for the single (uncoupled) neuromechanical module
*  single_oscillator_Eulerstep.m - gives an Euler's-method solution step for the single (uncoupled) neuromechanical module
* single_oscillator_PRC.m - computes the infinitessimal Phase Response Curves for the single (uncoupled) neuromechanical module
*  single_oscillator_adjoint_Eulerstep.m - gives an Euler's-method solution step for the adjoint equations to the single (uncoupled) neuromechanical module (needed to find the PRC)
* full_timetrace_to_phasediffs.m - computes phase-differences between each oscillator module in the full-state model over time from states
* full_timetrace_to_phasediffs.m - computes phase-differences between each oscillator module in the full-state model over time from states
* phasediffs_to_full_timetrace.m -  computes time-trace for the full state model from a vector of phase-locked phase differences between each oscillator module
* phases_to_init_cond.m - creates an initial condition for the full state model from a vector of initial phase differences between each oscillator module
* phase_locked_state_solve.m - solves for the phase-locked phase-differences of the full model using the phase model ODES (up to 4 modules)
 

### Figures/Plot Codes
* timetrace_to_curv_kym.m - creates a curvature kymograph from the time-trace of the full-state model
* blueblackred.m, bluewhitered.m, colorblind_colormap.mat, fireice.m, linspecer.m, shade.m - colormaps used for figures (not my own codes)

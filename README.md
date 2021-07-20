# MATLAB Code for two publications:

Carter L. Johnson (2020) "Neuromechanical Mechanisms of Locomotion in C. elegans" Ph.D. Dissertation, UC Davis. DOI: 10.13140/RG.2.2.24387.02082

Carter L. Johnson, Timothy J. Lewis, and Robert D. Guy (2021), "Neuromechanical Mechanisms of Gait Adaptation in C. elegans: Relative Roles of Neural and Mechanical Coupling", SIAM Journal on Applied Dynamical Systems. DOI: 10.1137/20M1346122

## Folders

* ./ - Main codes to run simulations, parameter sweeps, generate data and figures
* ./src - Source codes and functions that main codes run off of

 <details>
 <summary>How to Generate Figures from 2021 Paper</summary>
  
### How to Generate Figures from 2021 Paper

* Figure 3.1 - Run full_model_gammaloop.m
* Figure 3.2 - Run full_model_single_test_Run.m at given parameter values.
* Figures 3.3 and 3.4 - Run paperfig_parametrize_tmtf_loop_regimes.m to generate data.  Then run make_mukb_trend_colormap_figs.m to make fig.
* Figure 3.5 - Run full_model_eps_gamma_loop.m (commented out half generates data, uncommented half generates figs).
* Figure 4.1 - Run is_coupling_weak_single_module_vs_fully_coupled_module.m
* Figures 4.2, 4.3 and 4.4 - Run twobox_invest_loop.m
* Figure 4.5 - Run compare_gfunctions.m 
* Figure 4.6 - Run sixbox_osc_vs_full_gamma_loop.m
* Figure 4.7 - Run kymos_from_phasediffs_fig_sect4.m
* Figure C.1a - Run sixbox_osc_vs_full_gamma_loop_backwards_prop.m (commented out half generates data, uncommented half generates figs).
* Figures C.1b,c - Run twobox_invest_loop_backwardsprop.m
 </details>
 
  <details>
 <summary>How to Generate Figures from 2020 PhD Dissertation</summary>
 
### How to Generate Figures from 2020 PhD Dissertation

* Figure 2.2 -  Run full_model_gammaloop.m
* Figure 2.3 -  Run full_model_single_test_Run.m at given parameter values.
* Figures 2.4 and 2.5 - Run paperfig_parametrize_tmtf_loop_regimes.m to generate data.  Then run make_mukb_trend_colormap_figs.m to make fig.
* Figure 2.6 - Run full_model_eps_gamma_loop.m (commented out portions).
* Figure 3.1 - Run is_coupling_weak_single_module_vs_fully_coupled_module.m
* Figures 3.2, 3.3 and 3.4 - Run twobox_invest_loop.m
* Figure 3.5 - Run compare_gfunctions.m 
* Figure 3.6 - Run sixbox_osc_vs_full_gamma_loop.m
* Figure 3.7 - Run full_model_eps_gamma_loop_zero_muf.m
* Figures 4.1 and 4.2 - Bifurcation diagrams generated in [XPP/Auto](http://www.math.pitt.edu/~bard/bardware/tut/xppauto.html) (neuromechanical_module.ode), labels and colors added in Gimp.
* Figures 4.3 and 4.4 - Stability diagrams generated in [XPP/Auto](http://www.math.pitt.edu/~bard/bardware/tut/xppauto.html) (neuromechanical_module.ode), labels and colors added in Gimp.
* Figures 4.5 and 4.6 - Run oscillator_properties_contours.m
* Figure 4.7 - Stability diagrams generated in [XPP/Auto](http://www.math.pitt.edu/~bard/bardware/tut/xppauto.html) (neuromechanical_module_timescale_investigation.ode), labels and colors added in Gimp.
* Figures 4.8-4.13  - Run phase_response_properties.m
* Figures 4.14 and 4.15 - Run compare_PRCs_timescale_ordering.m
* Figure 5.2 - Phase plane created in [pplane](https://www.cs.unm.edu/~joel/dfield/), labels and colors added in Gimp.
* Figure 5.5 - Run find_alphas.m
* Figure 5.6 - Run period_amp_oscillations_from_1d_map.m
* Figures 6.1, 6.2 and 6.3 - Run diff_box_nos_phasemodels_vs_full_gamma_loop_fixed_ell.m to generate the data. Run diff_box_nos_phasemodels_vs_full_figures_fixed_ell.m to generate the figures.
* Figure 6.4 - Run diff_eps_20_box_phasemodels_vs_full_gamma_loop_fixed_ell.m
</details>

<details>
 <summary>Summary of src Codes</summary>
### Model Codes
* full_state_model_odes.m - function that generates ODEs for the forward locomotion system (dim neuromechanical modules)
  * full_state_model_odes_backwardsprop.m - same as above but with opposite-direction and signed proprioception
* coupled_oscillator_phase_difference_odes.m - function that generates phase-equation ODEs for the forward locomotion system (n oscillators)
   * coupled_oscillator_phase_difference_odes_backwardsprop.m - same as above but with opposite-direction and signed proprioception
* coup_currents_and_oscillator_coupling_fns.m / oscillator_coupling_fns.m - computes the coupling functions for the full coupled system from the single oscillator limit cycle and PRC
* single_oscillator_LC.m - computes the limit cycle oscillations for the single (uncoupled) neuromechanical module
*  single_oscillator_Eulerstep.m - gives an Euler's-method solution step for the single (uncoupled) neuromechanical module
* single_oscillator_PRC.m - computes the infinitessimal Phase Response Curves for the single (uncoupled) neuromechanical module
*  single_oscillator_adjoint_Eulerstep.m - gives an Euler's-method solution step for the adjoint equations to the single (uncoupled) neuromechanical module (needed to find the PRC)
* full_timetrace_to_phasediffs.m - computes phase-differences between each oscillator module in the full-state model over time from states
* full_timetrace_to_phasediffs.m - computes phase-differences between each oscillator module in the full-state model over time from states
* phasediffs_to_full_timetrace.m -  computes time-trace for the full state model from a vector of phase-locked phase differences between each oscillator module
* phases_to_init_cond.m - creates an initial condition for the full state model from a vector of initial phase differences between each oscillator module
* phase_locked_state_solve.m - solves for the phase-locked phase-differences of the full model using the phase model ODES (up to dim=4 modules)
* neuromechanical_module.ode - [XPP/Auto](http://www.math.pitt.edu/~bard/bardware/tut/xppauto.html) ode file for the single neuromechanical oscillator module 
 

### Figures/Plot Codes
* timetrace_to_curv_kym.m - creates a curvature kymograph from the time-trace of the full-state model
* blueblackred.m, bluewhitered.m, colorblind_colormap.mat, fireice.m, linspecer.m, shade.m - colormaps used for figures (not my own codes)
  </details>

# Stochastic model for site-specific amplification
Tools suite for characterizing seismic ground motion amplification 
between ground surface and reference points (depth or rock outcrop).
***************************************

This repository provides a comprehensive tools suite for the computation 
of a **Stochastic Model (SM)** designed to quantify ground motion 
amplification. The methodology accounts for uncertainties in subsurface 
properties by utilizing randomly perturbed 1D velocity models. The 
amplification is evaluated in the spectral domain and characterized by:
*   S/B Spectral Ratio (Amplification)
*   Energy Spectral Density Ratio (expressed in dB)
*   Envelope Delay (seconds)

**Applications:**
*   **Seismic Hazard Assessment:** Site-specific characterization for urban areas and critical infrastructure.
*   **Nuclear Safety:** Ground motion prediction for deep geological repositories (nuclear waste storage).
*   **Catastrophe Modeling:** Quantification of site effects and spectral amplification for risk assessment.
*   **Waveform Prediction:** Full-waveform broadband predictions at various depths.

1 METHODOLOGY
===================

The core of the toolset is based on the computation of 1D transfer functions. By incorporating 
stochastic perturbations of seismic velocities and layer thicknesses, the model provides a 
probabilistic view of site response, moving beyond simple deterministic estimates to full 
Uncertainty Quantification.

  Hallo, M., Bergamo, P., Fäh, D. (2022). Stochastic model to characterize 
high-frequency ground motion at depth validated by KiK-net vertical array data,
Bulletin of the Seismological Society of America, 112 (4), 1997–2017. [https://doi.org/10.1785/0120220038](https://doi.org/10.1785/0120220038)

2 TECHNICAL IMPLEMENTATION
===================

Cross-Platform (Windows, Linux), Automatic exports

3 PACKAGE CONTENT
===================

  1. `SM.m` - Subroutine for computation of the Stochastic Model (SM)
  2. `RUNSM.m` - Example of wrapper code to run the SM subroutine

4 REQUIREMENTS
===================

  MATLAB: Version R2018b, Codes do not require any additional Matlab Toolboxes.

5 USAGE
===================

  1. Open MATLAB
  2. Run the main scripts: `RUNSM.m`

6 COPYRIGHT
===================

Copyright (C) 2020-2022 Swiss Seismological Service, ETH Zurich

This program is published under the GNU General Public License (GNU GPL).

This program is free software: you can modify it and/or redistribute it
or any derivative version under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version.

This code is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
and don't remove their names from the code.

You should have received copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

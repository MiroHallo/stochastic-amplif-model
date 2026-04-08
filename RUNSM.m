%%
% It computes the Stochastic Model (SM) of the ground motion amplification
% for a given local 1D velocity model (1D resonance effects and impedance)
%
% Hallo, M., Bergamo, P., Fäh, D. (2022). Stochastic model to characterize
%      high-frequency ground motion at depth validated by KiK-net vertical
%      array data, Bulletin of the Seismological Society of America.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code author: Miroslav Hallo
% ETH Zurich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Revision 9/2021: The first version of the wrapper function
% Revision 5/2022: Adjusted for distribution
% Tested in Matlab R2018b, R2021a
%
% Copyright (C) 2021,2022  Swiss Seismological Service, ETH Zurich
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT:
close all;
clear all;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% Velocity model
Depth0 = [0 4 8 10 24 52]; % Top depth of layers [m]
Vs0 = [131 183 183 386 1076 1469]; % Shear wave velocity [m/s]
Rho0 = [1700 2000 2000 2200 2200 2200]; % Density [kg/m3]
Qs0 = [13 18 18 39 108 147]; % Shear wave quality factor Q [-]
Fox0 = 1./(2*Qs0); % Layer damping ratio

% Reference depth [m]
% (depth of the borehole receiver or outcrop base-rock)
RefDepth = 110;

% Number of perturbed models for construction of the Stochastic Model
% (recommended: 100-10000)
nM = 1000;

% Percentage of velocity model perturbation 1sigma [%]
% (shear modulus, density, damping, interface depth)
% Set 1sigma of density into -1 (e.g. [15,-1,15,15]) for correlated
% densities from empirical relation by Nagashima & Kawase (2021)
pert1s = [15,15,15,15];

% 1sigma [deg] of the incident angle distribution
% (centered on the vertical incidence angle, recommended: 5-25)
inc1s = 15;

% Frequencies [Hz] of the output spectral curve
% (min, max, samples)
fMinMaxN = [0.1,50,1000];

% Compute transfer functions for Borehole or Rock Outcrop
% (Outcrop=0, Borehole=1, Incidence=2)
BorO = 1;

% Plot figures?
% (No = 0, Models = 1, All = 2)
PlotFigs = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the Stochastic Model
[ESD_mean,ESD_sigma,ED_mean_log,ED_sigma_log] = SM(Depth0,Vs0,Rho0,Fox0,...
    RefDepth,nM,pert1s,inc1s,fMinMaxN,BorO,PlotFigs);

% Resultant Spectral ratio [-]
SBratio_mean = sqrt(10.^(ESD_mean./10));
SBratio_upper = sqrt(10.^((ESD_mean+ESD_sigma)./10));
SBratio_lower = sqrt(10.^((ESD_mean-ESD_sigma)./10));

% Resultant Energy Spectral Density ratio [dB]
ESD_mean;
ESD_upper = ESD_mean+ESD_sigma;
ESD_lower = ESD_mean-ESD_sigma;

% Resultant Envelope Delay [s]
ED_mean = 10.^ED_mean_log;
ED_upper = 10.^(ED_mean_log+ED_sigma_log);
ED_lower = 10.^(ED_mean_log-ED_sigma_log);

% Save the Stochastic Model into ASCII file
f = linspace(fMinMaxN(1),fMinMaxN(2),fMinMaxN(3));
fid = fopen('Output.dat','w');
fprintf(fid,'%s\r\n','# Stochastic Model following Hallo et al. (2022)');
fprintf(fid,'%s\r\n','# Frequency | S/B ratio | S/B ESD rat. | S-B delay');
fprintf(fid,'%s\r\n','# Hz        | -         | dB           | sec');
for i=1:length(f)
    fprintf(fid,'%11.5f %11.5f %11.5f %11.5f\r\n',...
        [f(i),SBratio_mean(i),ESD_mean(i),ED_mean(i)]);
end
fclose(fid);


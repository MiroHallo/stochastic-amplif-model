function [ESD_mean, ESD_sigma, ED_mean_log, ED_sigma_log] = ...
    sm(Depth0, Vs0, Rho0, Fox0, RefDepth, nM, pert1s, inc1s, fMinMaxN, BorO, PlotFigs)
% SM Stochastic Model to characterize site-specific amplification
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Compute the Stochastic Model (SM) to characterize seismic ground motion 
% amplification between ground surface and reference points (depth or rock 
% outcrop) including stochastic perturbations of local 1D velocity models.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% ETH Zürich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Revision 2020/03: The first version of the function
% Revision 2021/09: Updated after experiments with real data
% Revision 2022/05: Adjusted for distribution
% Revision 2023/02: New features for the TransferC purposes
% Tested in Matlab R2018b, R2021a, R2022a, R2025b
% Method:
% Hallo, M., Bergamo, P., Fäh, D. (2022). Stochastic model to characterize
%      high-frequency ground motion at depth validated by KiK-net vertical
%      array data, Bulletin of the Seismological Society of America, 112(4),
%      1997–2017. https://doi.org/10.1785/0120220038
%
% Copyright (C) 2020-2023  Swiss Seismological Service, ETH Zurich
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
%
% INPUT:
% Depth0 - Top depth of layers [m] (vector)
% Vs0 - Shear wave velocity of layers [m/s] (vector)
% Rho0 - Density of layers [kg/m3] (vector)
% Fox0 - Damping ratio of layers (vector)
% RefDepth - Reference depth [m]
% nM - Number of perturbed models
% pert1s - Perturbation [%] (shear modulus, density, damping, interfaces)
%          It can be vector (joint for all layers) or matrix (individual)
%          Set pert1s(1,1) negative for a direct perturbation of velocity
%          Set pert1s(1,2) negative for densities from Nagashima & Kawase (2021)
% inc1s - Incident angle distribution [deg] (1sigma centered on vertical)
% fMinMaxN - Frequencies [Hz] of the output curve (min, max, samples)
% BorO - Borehole or Rock Outcrop? (Outcrop=0, Borehole=1, Incidence=2)
% PlotFigs - Plot figures? (No = 0, Models = 1, All = 2)
%
% OUTPUT:
% ESD_mean - S/B Energy Spectral Density ratio [dB]
% ESD_sigma - S/B Energy Spectral Density error (1sigma)
% ED_mean_log - S-B Envelope Delay [log10(s)]
% ED_sigma_log - S-B Envelope Delay error (1sigma)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ESD_mean = NaN;
ESD_sigma = NaN;
ED_mean_log = NaN;
ED_sigma_log = NaN;

% Create row-vectors
Depth0 = Depth0(:)';
Vs0 = Vs0(:)';
Rho0 = Rho0(:)';
Fox0 = Fox0(:)';

% Number of layers
nL = length(Depth0);
if (length(Vs0) ~= nL) || (length(Rho0) ~= nL) || (length(Fox0) ~= nL)
    disp('ERROR: Input model vectors are not of the same length')
    return
end

% Reference depth has to be below all layers
IR = find(Depth0 >= RefDepth,1,'first');
if ~isempty(IR)
    if IR==1
        disp('ERROR: Reference depth is above the shallowest interface')
        return
    else
        % Remove layers below the reference depth
        Depth0 = Depth0(1:IR-1);
        Vs0 = Vs0(1:IR-1);
        Rho0 = Rho0(1:IR-1);
        Fox0 = Fox0(1:IR-1);
        if min(size(pert1s)) > 1
            pert1s = pert1s(1:IR-1,:);
        end
        nL = length(Depth0);
    end
end

% Number of perturbed models
if nM<3
    disp('ERROR: Number of perturbed models must be >2')
    return
end

inc1s = abs(inc1s);
if inc1s>30
    disp('ERROR: 1sigma of the incident angle has to be <=30 [deg]')
    return
end

% Prepare frequencies
f = linspace(fMinMaxN(1), fMinMaxN(2), fMinMaxN(3));

% Prepare switches
if pert1s(1,1) < 0
    nkvelo = 1;
    pert1s(:,1) = abs(pert1s(:,1));
else
    nkvelo = 0;
end

if pert1s(1,2) < 0
    nkdens = 1;
    pert1s(:,2) = pert1s(:,1);
else
    nkdens = 0;
end

perDim = size(pert1s);
if min(perDim) == 1
    pert1s = repmat(pert1s,nL,1);
end

% convert percentage into probability
pert1sP = abs(pert1s) ./ 100;


%% -------------------------------------------------------------------
% Step 1 - Perturb velocity models

% Allocate
Depth = zeros(nM,nL);
Vs = zeros(nM,nL);
Rho = zeros(nM,nL);
Fox = zeros(nM,nL);

% The first model is the initial one
Depth(1,1:nL) = Depth0;
Vs(1,1:nL) = Vs0;
Rho(1,1:nL) = Rho0;
Fox(1,1:nL) = Fox0;

% Get shear modulus G0
G0 = (Vs0.^2).*Rho0;

% Create perturbed velocity models
for mod = 2:nM
    % Perturb density
    Rho(mod,:) = Rho0 .* exp(randn(1,nL).*pert1sP(1:nL,2)');
    % Perturb damping
    Fox(mod,:) = Fox0 .* exp(randn(1,nL).*pert1sP(1:nL,3)');
    if nkvelo == 1
        % Perturb S-wave velocity
        Vs(mod,:) = Vs0 .* exp(randn(1,nL).*pert1sP(1:nL,1)');
    else
        % Perturb shear modulus
        G = G0 .* exp(randn(1,nL).*pert1sP(1:nL,1)');
        Vs(mod,:) = sqrt(G./Rho(mod,:));
    end
    % Perturb depths of interfaces
    OKtest = 0;
    count = 0;
    pertTmp = pert1sP(1:nL,4)';
    while OKtest ~= 1
        count = count + 1;
        Z = sqrt((Depth0.^2).*exp(randn(1,nL).*pertTmp));
        if ~sum(sign(diff(Z))<0) && ~sum(Z<0)
            OKtest = 1;
        end
        if count > 100
            pertTmp = pertTmp./2;
            count = 0;
        end
    end
    Depth(mod,:) = Z;
end

% Densities from empirical relation by Nagashima & Kawase (2021)
if nkdens == 1
    for mod = 1:nM
        Rho(mod,:) = 1000.*(1+0.187.*(Vs(mod,:).^0.263));
    end
end


%% -------------------------------------------------------------------
% Step 2 - Ensemble of Transfer Functions

% Incident angle perturbations
if inc1s==0
    th_all = zeros(nM,1);
else
    th_dist = truncate(makedist('Normal',0,inc1s),-60,60);
    th_all = abs(random(th_dist,nM,1));
    th_all(1,1) = 0;
end

% Compute SH-wave Transfer Functions
Nf = length(f);
df = f(2)-f(1);
h_all = zeros(Nf,nM);
ed_all = zeros(Nf,nM);
for mod = 1:nM
    [h1,h2,h3] = respSH(Depth(mod,:),Vs(mod,:),Rho(mod,:),Fox(mod,:),f,th_all(mod));

    % Chose Transfer Function type (Outcrop=0, Borehole=1, Incidence=2)
    if BorO == 1
        h0 = h2;
    elseif BorO == 2
        h0 = h1;
    else
        h0 = h3;
    end
    % Final Transfer Function and ED
    Y = conj(h0);
    dfi = [0, diff(unwrap(angle(Y))) / df];
    h_all(:,mod) = h0;
    ed_all(:,mod) = dfi/(2*pi);
end


%% -------------------------------------------------------------------
% Step 3 - Stochastic Model

% Compute mean spectral curves
ESD_all = abs(h_all).^2;
ESD_all = 10.*log10(ESD_all);
ESD_mean = mean(ESD_all,2); % Energy spectral density (ESD)
ED_mean_log = mean(log10(ed_all),2); % Envelope delay (ED)

% Compute variances of ESD and ED
ESD_var = zeros(Nf,1);
ED_var_log = zeros(Nf,1);
for j = 1:Nf
    for n = 1:nM
        ESD_var(j) = ESD_var(j) + (ESD_all(j,n) - ESD_mean(j))^2;
        ED_var_log(j) = ED_var_log(j) + (log10(ed_all(j,n)) - ED_mean_log(j))^2;
    end
end
ESD_var = ESD_var./nM;
ED_var_log = ED_var_log./nM;

% Compute 1sigma of ESD and ED
ESD_sigma = sqrt(ESD_var);
ED_sigma_log = sqrt(ED_var_log);
ED_sigma_log(isnan(ED_sigma_log)) = 0;


%% -------------------------------------------------------------------
% Plot figure 1

if PlotFigs < 1
    return
end

freqlim = [min(f),max(f)];

% Prepare ticks
if max(f)>10
    fTick = [0.01 0.1 1 10 max(f)];
    fTicklabel = {'0.01' '0.1' '1' '10' num2str(round(max(f)))};
else
    fTick = freqlim;
    fTicklabel = {num2str(round(min(f))) num2str(round(max(f)))};
end

% Prepare text
if BorO == 1
    rr = 'B';
elseif BorO == 2
    rr = 'I';
else
    rr = 'R';
end

% Threshold for models to plot
nM2 = min(nM,500);

% Prepare layered models for plots
lay2plot = zeros(4,nM2,2*nL);
for mod = 1:nM2
    lay2plot(1,mod,1:2:2*nL) = Depth(mod,1:end);
    lay2plot(1,mod,2:2:(2*nL)-2) = Depth(mod,2:end);
    lay2plot(1,mod,2*nL) = RefDepth;
    lay2plot(2,mod,1:2:2*nL) = Vs(mod,1:end);
    lay2plot(2,mod,2:2:2*nL) = Vs(mod,1:end);
    lay2plot(3,mod,1:2:2*nL) = Rho(mod,1:end);
    lay2plot(3,mod,2:2:2*nL) = Rho(mod,1:end);
    lay2plot(4,mod,1:2:2*nL) = Fox(mod,1:end);
    lay2plot(4,mod,2:2:2*nL) = Fox(mod,1:end);
end

% Perturbed velocity models
figure('Color','w','Units','normalized','OuterPosition',[0.2,0.1,0.6,0.8]);
cc = colormap(['copper(',num2str(nM2),')']);
randcolol = max(1,(1:nM2)' - randi(floor(nM2/3),nM2,1));
cc = cc(randcolol,:);

subplot(1,4,1:2)
hln=[];
for mod = 2:nM2
    hln(2)=plot(squeeze(lay2plot(2,mod,:)),squeeze(lay2plot(1,mod,:)),'Color',cc(mod,:));
    hold on;
end
hln(1)=plot(squeeze(lay2plot(2,1,:)),squeeze(lay2plot(1,1,:)),'k--','LineWidth',1.5); 
legend(hln,'InModel','Ensemble','Location','northeast')
hold off;
set(gca,'Layer','top');
set(gca,'YDir','reverse'); set(gca,'YLim',[0 RefDepth]);
xlabel('S-wave velocity (m/s)');
ylabel('Depth (m)');
set(gca,'FontSize',10); box on;

subplot(1,4,3)
for mod = 2:nM2
    plot(squeeze(lay2plot(3,mod,:)),squeeze(lay2plot(1,mod,:)),'Color',cc(mod,:));
    hold on;
end
plot(squeeze(lay2plot(3,1,:)),squeeze(lay2plot(1,1,:)),'k--','LineWidth',1.5);
hold off;
set(gca,'Layer','top');
set(gca,'YDir','reverse');
set(gca,'YLim',[0 RefDepth]);
xlabel('Density (kg/m^3)');
ylabel('Depth (m)');
set(gca,'FontSize',10); box on;

subplot(1,4,4)
for mod = 2:nM2
    plot(squeeze(lay2plot(4,mod,:)),squeeze(lay2plot(1,mod,:)),'Color',cc(mod,:));
    hold on;
end
plot(squeeze(lay2plot(4,1,:)),squeeze(lay2plot(1,1,:)),'k--','LineWidth',1.5);
hold off;
set(gca,'Layer','top');
set(gca,'YDir','reverse');
set(gca,'YLim',[0 RefDepth]);
xlabel('Attenuation \xi');
ylabel('Depth (m)');
set(gca,'FontSize',10); box on;


%% -------------------------------------------------------------------
% Plot figure 2

if PlotFigs < 2
    return
end

% Ensemble of Transfer Functions
figure('Color','w','Units','normalized','OuterPosition',[0.2,0.1,0.6,0.8]);

subplot(3,2,1)
hln=[];
for mod = 2:nM2
    hln(3) = semilogx(f,abs(h_all(:,mod)),'Color',cc(mod,:)); hold on;
end
hln(2) = semilogx(f,abs(h_all(:,1)),'k--','Linewidth',1.5);
hln(1) = semilogx(f,sqrt(10.^(ESD_mean./10)),'Color','b','Linewidth',2);
hold off;
legend(hln,'SM','InModel','Ensemble','Location','northwest');
set(gca,'Xtick',fTick);
set(gca,'XtickLabel',fTicklabel);
set(gca,'Xlim',freqlim);
xtickangle(0);
set(gca,'Layer','top');
xlabel('Frequency (Hz)')
ylabel(['S/',rr,' spectral ratio'])
set(gca,'FontSize',10); box on;

subplot(3,2,3)
for mod = 2:nM2
    semilogx(f,ESD_all(:,mod),'Color',cc(mod,:)); hold on;
end
semilogx(f,ESD_all(:,1),'k--','Linewidth',1.5);
semilogx(f,ESD_mean,'Color','b','Linewidth',2);
hold off;
set(gca,'Xtick',fTick);
set(gca,'XtickLabel',fTicklabel);
set(gca,'Xlim',freqlim);
xtickangle(0);
set(gca,'Layer','top');
xlabel('Frequency (Hz)')
ylabel(['S/',rr,' ESD ratio (dB)'])
set(gca,'FontSize',10); box on;

subplot(3,2,5)
for mod = 2:nM2
    loglog(f,ed_all(:,mod),'Color',cc(mod,:)); hold on;
end
loglog(f,ed_all(:,1),'k--','Linewidth',1.5);
loglog(f,10.^ED_mean_log,'Color','b','Linewidth',2);
hold off;
set(gca,'Xtick',fTick);
set(gca,'XtickLabel',fTicklabel);
set(gca,'Xlim',freqlim);
xtickangle(0);
set(gca,'Layer','top');
xlabel('Frequency (Hz)')
ylabel(['S-',rr,' envelope delay (s)'])
set(gca,'FontSize',10); box on;

% -------------------------------------------------------------------
% Stochastic Model

subplot(3,2,2)
ESD_tmp = sqrt(10.^(ESD_mean./10));
hln=[];
hln(1) = semilogx(f,ESD_tmp,'Color','b','Linewidth',2); hold on;
ESD_tmp = sqrt(10.^((ESD_mean+ESD_sigma)./10));
hln(2) = semilogx(f,ESD_tmp,'--','Color','b','Linewidth',1);
ESD_tmp = sqrt(10.^((ESD_mean-ESD_sigma)./10));
semilogx(f,ESD_tmp,'--','Color','b','Linewidth',1);
ESD_tmp = sqrt(10.^((ESD_mean+2*ESD_sigma)./10));
hln(3) = semilogx(f,ESD_tmp,':','Color','b','Linewidth',0.8);
ESD_tmp = sqrt(10.^((ESD_mean-2*ESD_sigma)./10));
semilogx(f,ESD_tmp,':','Color','b','Linewidth',0.8);
hold off;
legend(hln,'SM','Error 1\sigma','Error 2\sigma','Location','northwest');
set(gca,'Xtick',fTick);
set(gca,'XtickLabel',fTicklabel);
set(gca,'Xlim',freqlim);
xtickangle(0); set(gca,'Layer','top');
xlabel('Frequency (Hz)')
ylabel(['S/',rr,' spectral ratio'])
set(gca,'FontSize',10); box on;

subplot(3,2,4)
ESD_tmp = ESD_mean;
semilogx(f,ESD_tmp,'Color','b','Linewidth',2); hold on;
ESD_tmp = ESD_mean+ESD_sigma;
semilogx(f,ESD_tmp,'--','Color','b','Linewidth',1);
ESD_tmp = ESD_mean-ESD_sigma;
semilogx(f,ESD_tmp,'--','Color','b','Linewidth',1);
ESD_tmp = ESD_mean+2*ESD_sigma;
semilogx(f,ESD_tmp,':','Color','b','Linewidth',0.8);
ESD_tmp = ESD_mean-2*ESD_sigma;
semilogx(f,ESD_tmp,':','Color','b','Linewidth',0.8);
hold off;
set(gca,'Xtick',fTick);
set(gca,'XtickLabel',fTicklabel);
set(gca,'Xlim',freqlim);
xtickangle(0); set(gca,'Layer','top');
xlabel('Frequency (Hz)')
ylabel(['S/',rr,' ESD ratio (dB)'])
set(gca,'FontSize',10); box on;

subplot(3,2,6)
ED_tmp = 10.^ED_mean_log;
loglog(f,ED_tmp,'Color','b','Linewidth',2); hold on;
ED_tmp = 10.^(ED_mean_log+ED_sigma_log);
loglog(f,ED_tmp,'--','Color','b','Linewidth',1);
ED_tmp = 10.^(ED_mean_log-ED_sigma_log);
loglog(f,ED_tmp,'--','Color','b','Linewidth',1);
ED_tmp = 10.^(ED_mean_log+2*ED_sigma_log);
loglog(f,ED_tmp,':','Color','b','Linewidth',0.8);
ED_tmp = 10.^(ED_mean_log-2*ED_sigma_log);
loglog(f,ED_tmp,':','Color','b','Linewidth',0.8);
hold off;
set(gca,'Xtick',fTick);
set(gca,'XtickLabel',fTicklabel);
set(gca,'Xlim',freqlim);
xtickangle(0); set(gca,'Layer','top');
xlabel('Frequency (Hz)')
ylabel(['S-',rr,' envelope delay (s)'])
set(gca,'FontSize',10); box on;

end


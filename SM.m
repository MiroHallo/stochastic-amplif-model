function [ESDmean,ESDsigma,EDlmean,EDlsigma] = SM(Depth0,Vs0,Rho0,Fox0,...
    RefDepth,nM,pert1s,inc1s,fMinMaxN,BorO,PlotFigs)
%
% Compute the Stochastic Model (SM):
% Hallo, M., Bergamo, P., Fäh, D. (2022). Stochastic model to characterize
%      high-frequency ground motion at depth validated by KiK-net vertical
%      array data, Bulletin of the Seismological Society of America.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code author: Miroslav Hallo
% ETH Zurich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Revision 3/2020: The first version of the function
% Revision 9/2021: Updated after experiments with real data
% Revision 5/2022: Adjusted for distribution
% Revision 2/2023: New features for the TransferC purposes
% Tested in Matlab R2018b, R2021a, R2022a
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
% inc1s - Incident angle distribution [deg] (1sigma centered on vertical)
% fMinMaxN - Frequencies [Hz] of the output curve (min, max, samples)
% BorO - Borehole or Rock Outcrop? (Outcrop=0, Borehole=1, Incidence=2)
% PlotFigs - Plot figures? (No = 0, Models = 1, All = 2)
%
% OUTPUT:
% ESDmean - S/B Energy Spectral Density ratio [dB]
% ESDsigma - S/B Energy Spectral Density error (1sigma)
% EDlmean - S-B Envelope Delay [log10(s)]
% EDlsigma - S-B Envelope Delay error (1sigma)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ESDmean=NaN; ESDsigma=NaN; EDlmean=NaN; EDlsigma=NaN;
Depth0=Depth0(:)'; Vs0=Vs0(:)'; Rho0=Rho0(:)'; Fox0=Fox0(:)';
nL=length(Depth0);
if (length(Vs0)~=nL)||(length(Rho0)~=nL)||(length(Fox0)~=nL)
    disp('ERROR: Input model vectors are not of the same length')
    return
end
IR=find(Depth0>=RefDepth,1,'first');
if ~isempty(IR)
    if IR==1
        disp('ERROR: Reference depth is above the shallowest interface')
        return
    else
        Depth0=Depth0(1:IR-1);
        Vs0=Vs0(1:IR-1);
        Rho0=Rho0(1:IR-1);
        Fox0=Fox0(1:IR-1);
        if min(size(pert1s))>1
            pert1s=pert1s(1:IR-1,:);
        end
        nL=length(Depth0);
    end
end
if nM<3
    disp('ERROR: Number of perturbed models must be >2')
    return
end
inc1s=abs(inc1s);
if inc1s>30
    disp('ERROR: 1sigma of the incident angle dist. has to be <=30 [deg]')
    return
end
f=linspace(fMinMaxN(1),fMinMaxN(2),fMinMaxN(3));
if pert1s(1,1)<0, nkvelo=1; pert1s(:,1)=abs(pert1s(:,1));
else, nkvelo=0; end
if pert1s(1,2)<0, nkdens=1; pert1s(:,2)=pert1s(:,1);
else, nkdens=0; end
perDim=size(pert1s);
if min(perDim)==1
    pert1s=repmat(pert1s,nL,1);
end
pert1sP=abs(pert1s)./100;

% Compute the Stochastic Model
Depth=zeros(nM,nL); Depth(1,1:nL)=Depth0;
Vs=zeros(nM,nL); Vs(1,1:nL)=Vs0;
Rho=zeros(nM,nL); Rho(1,1:nL)=Rho0;
Fox=zeros(nM,nL); Fox(1,1:nL)=Fox0;
G0=(Vs0.^2).*Rho0;
for mod = 2:nM
    Rho(mod,:)=Rho0.*exp(randn(1,nL).*pert1sP(1:nL,2)');
    Fox(mod,:)=Fox0.*exp(randn(1,nL).*pert1sP(1:nL,3)');
    if nkvelo==1
        Vs(mod,:)=Vs0.*exp(randn(1,nL).*pert1sP(1:nL,1)');
    else
        G=G0.*exp(randn(1,nL).*pert1sP(1:nL,1)');
        Vs(mod,:)=sqrt(G./Rho(mod,:));
    end
    Okt=0;
    count=0;
    pertTmp=pert1sP(1:nL,4)';
    while Okt ~= 1
        count=count+1;
        Z=sqrt((Depth0.^2).*exp(randn(1,nL).*pertTmp));
        if ~sum(sign(diff(Z))<0) && ~sum(Z<0), Okt = 1; end
        if count>100, pertTmp=pertTmp./2; count = 0; end
    end
    Depth(mod,:)=Z;
end
if nkdens==1
    for mod = 1:nM
        Rho(mod,:)=1000.*(1+0.187.*(Vs(mod,:).^0.263));
    end
end
if inc1s==0, uda=zeros(nM,1);
else
    udad=truncate(makedist('Normal',0,inc1s),-60,60);
    uda=abs(random(udad,nM,1));
    uda(1,1)=0;
end
Nf=length(f); df=f(2)-f(1);
h1=zeros(Nf,nM);
ed1=zeros(Nf,nM);
for mod = 1:nM
    Thi=zeros(1,nL);
    Thi(1:end-1)=Depth(mod,2:end)-Depth(mod,1:end-1); Thi(end)=Inf;
    Vsi=Vs(mod,:).*(1+1i*Fox(mod,:));
    Slh=sind(uda(mod))/Vs(mod,end);
    Slv=sqrt((1./(Vsi.^2))-(Slh^2));
    Vsp=(Vsi.^2).*Slv;
    Alp=(Rho(mod,1:end-1).*Vsp(1:end-1))./(Rho(mod,2:end).*Vsp(2:end));
    A=complex(ones(nL,Nf)); B=complex(ones(nL,Nf));
    for m = 1:nL-1
        Ksh=2*pi*f*Thi(m)*Slv(m);
        A(m+1,1:Nf)=0.5*A(m,1:Nf)*(1+Alp(m)).*exp(1i*Ksh)+ ...
            0.5*B(m,1:Nf)*(1-Alp(m)).*exp(-1i*Ksh);
        B(m+1,1:Nf)=0.5*A(m,1:Nf)*(1-Alp(m)).*exp(1i*Ksh)+ ...
            0.5*B(m,1:Nf)*(1+Alp(m)).*exp(-1i*Ksh);
    end
    if BorO==1
        h2(1:Nf)=(A(1,:)+B(1,:))./(A(nL,:)+B(nL,:));
    elseif BorO==2
        h2(1:Nf)=2./A(nL,:);
    else
        h2(1:Nf)=1./A(nL,:);
    end
    Y=conj(h2);
    dfi=[0,diff(unwrap(angle(Y)))/df];
    h1(:,mod)=h2;
    ed1(:,mod)=dfi/(2*pi);
end
h3=abs(h1).^2;
h3=10.*log10(h3);
ESDmean=mean(h3,2);
EDlmean=mean(log10(ed1),2);
h3v=zeros(Nf,1);
ed1v=zeros(Nf,1);
for j=1:Nf
    for n = 1:nM
        h3v(j)=h3v(j)+(h3(j,n)-ESDmean(j))^2;
        ed1v(j)=ed1v(j)+(log10(ed1(j,n))-EDlmean(j))^2;
    end
end
h3v=h3v./nM;
ed1v=ed1v./nM;
ESDsigma=sqrt(h3v);
EDlsigma=sqrt(ed1v);
EDlsigma(isnan(EDlsigma))=0;

% Prepare figures
if PlotFigs < 1
    return
end
freqlim=[min(f),max(f)];
if max(f)>10
    fTi=[0.01 0.1 1 10 max(f)];
    fTl={'0.01' '0.1' '1' '10' num2str(round(max(f)))};
else
    fTi=freqlim;
    fTl={num2str(round(min(f))) num2str(round(max(f)))};
end
if BorO==1
    rr = 'B';
elseif BorO==2
    rr = 'In';
else
    rr = 'O';
end
lm=zeros(4,min(nM,500),2*nL);
for mod = 1:min(nM,500)
    lm(1,mod,1:2:2*nL)=Depth(mod,1:end);
    lm(1,mod,2:2:(2*nL)-2)=Depth(mod,2:end); lm(1,mod,2*nL)=RefDepth;
    lm(2,mod,1:2:2*nL)=Vs(mod,1:end); lm(2,mod,2:2:2*nL)=Vs(mod,1:end);
    lm(3,mod,1:2:2*nL)=Rho(mod,1:end); lm(3,mod,2:2:2*nL)=Rho(mod,1:end);
    lm(4,mod,1:2:2*nL)=Fox(mod,1:end); lm(4,mod,2:2:2*nL)=Fox(mod,1:end);
end

% Plot perturbed velocity models
figure('Color','w','Units','normalized','OuterPosition',[0.2,0.1,0.6,0.8]);
cc=colormap(['copper(',num2str(min(nM,500)),')']);
ranc=max(1,(1:min(nM,500))'-randi(floor(min(nM,500)/3),min(nM,500),1));
cc=cc(ranc,:); u=[];
subplot(1,4,1:2)
for mod = 2:min(nM,500)
    u(2)=plot(squeeze(lm(2,mod,:)),squeeze(lm(1,mod,:)),'Color',cc(mod,:));
    hold on;
end
u(1)=plot(squeeze(lm(2,1,:)),squeeze(lm(1,1,:)),'k--','LineWidth',2);
legend(u,'InModel','Ensemble','Location','northeast','Interpreter','tex')
hold off; set(gca,'Layer','top');
set(gca,'YDir','reverse'); set(gca,'YLim',[0 RefDepth]);
xlabel('S-wave velocity (m/s)','Interpreter','tex');
ylabel('Depth (m)','Interpreter','tex');
set(gca,'FontSize',10); box on;
subplot(1,4,3)
for mod = 2:min(nM,500)
    plot(squeeze(lm(3,mod,:)),squeeze(lm(1,mod,:)),'Color',cc(mod,:));
    hold on;
end
plot(squeeze(lm(3,1,:)),squeeze(lm(1,1,:)),'k--','LineWidth',2);
hold off; set(gca,'Layer','top');
set(gca,'YDir','reverse'); set(gca,'YLim',[0 RefDepth]);
xlabel('Density (kg/m^3)','Interpreter','tex');
ylabel('Depth (m)','Interpreter','tex');
set(gca,'FontSize',10); box on;
subplot(1,4,4)
for mod = 2:min(nM,500)
    plot(squeeze(lm(4,mod,:)),squeeze(lm(1,mod,:)),'Color',cc(mod,:));
    hold on;
end
plot(squeeze(lm(4,1,:)),squeeze(lm(1,1,:)),'k--','LineWidth',2);
hold off; set(gca,'Layer','top');
set(gca,'YDir','reverse'); set(gca,'YLim',[0 RefDepth]);
xlabel('Damping \xi','Interpreter','tex');
ylabel('Depth (m)','Interpreter','tex');
set(gca,'FontSize',10); box on;

% Plot the Stochastic Model
if PlotFigs < 2
    return
end
figure('Color','w','Units','normalized','OuterPosition',[0.2,0.1,0.6,0.8]);
subplot(3,1,1)
ESDtmp=sqrt(10.^(ESDmean./10)); u=[];
u(1)=semilogx(f,ESDtmp,'Color','b','Linewidth',2); hold on;
ESDtmp=sqrt(10.^((ESDmean+ESDsigma)./10));
u(2)=semilogx(f,ESDtmp,':','Color','b','Linewidth',1);
ESDtmp=sqrt(10.^((ESDmean-ESDsigma)./10));
semilogx(f,ESDtmp,':','Color','b','Linewidth',1); hold off;
legend(u,'SM','Error 1\sigma','Location','northwest','Interpreter','tex');
set(gca,'Xtick',fTi); set(gca,'XtickLabel',fTl); set(gca,'Xlim',freqlim);
xtickangle(0); set(gca,'Layer','top');
xlabel('Frequency (Hz)','Interpreter','tex')
ylabel(['S/',rr,' spectral ratio'],'Interpreter','tex')
title('Stochastic Model','Interpreter','tex')
set(gca,'FontSize',10); box on;
subplot(3,1,2)
semilogx(f,ESDmean,'Color','b','Linewidth',2); hold on;
semilogx(f,(ESDmean+ESDsigma),':','Color','b','Linewidth',1);
semilogx(f,(ESDmean-ESDsigma),':','Color','b','Linewidth',1); hold off;
set(gca,'Xtick',fTi); set(gca,'XtickLabel',fTl); set(gca,'Xlim',freqlim);
xtickangle(0); set(gca,'Layer','top');
xlabel('Frequency (Hz)','Interpreter','tex')
ylabel(['S/',rr,' ESD ratio (dB)'],'Interpreter','tex')
set(gca,'FontSize',10); box on;
subplot(3,1,3)
EDtmp=10.^EDlmean;
loglog(f,EDtmp,'Color','b','Linewidth',2); hold on;
EDtmp=10.^(EDlmean+EDlsigma);
loglog(f,EDtmp,':','Color','b','Linewidth',1);
EDtmp=10.^(EDlmean-EDlsigma);
loglog(f,EDtmp,':','Color','b','Linewidth',1); hold off;
set(gca,'Xtick',fTi); set(gca,'XtickLabel',fTl); set(gca,'Xlim',freqlim);
xtickangle(0); set(gca,'Layer','top');
xlabel('Frequency (Hz)','Interpreter','tex')
ylabel(['S-',rr,' envelope delay (s)'],'Interpreter','tex')
set(gca,'FontSize',10); box on;

end



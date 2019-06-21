% Script to verify that the SS equations give the expected free pool
% formulae when no MT is present. These results are not included in the
% publication

tissuepars = init_tissue();

%%% Remove MT
tissuepars.semi.M0 = 0;
tissuepars.k = 0;


%%% Set up sequence

%%% generate some pulses
flips = deg2rad(1:2:50);
nf = length(flips);
tau = 2e-3;
TR = 5e-3;
b1_rms = 4;
delta = 7e3;
dt = 10e-6;

pulses = {};
b1sqrd = {};
b1pulse= {};
for ii=1:nf
    % 3 bands
    [pulses{ii,1},b1sqrd{ii,1},b1pulse{ii,1}] = gen_MB_pulse(flips(ii),tau,TR,b1_rms,delta,'3','alpha',3,'dt',dt); 
    % 2 bands
    [pulses{ii,2},b1sqrd{ii,2},b1pulse{ii,2}] = gen_MB_pulse(flips(ii),tau,TR,b1_rms,delta,'2+','alpha',3,'dt',dt);
end


%% bSSFP steady state

%%% Extra pars needed for SSFP
dphi=0;

% loop over the different flips
Mssfp = zeros(nf,3);
Delta_Hz = [-delta 0 delta];
for ii=1:nf
    %%% Symmetric (3-band) case
    Mssfp(ii,1) = ssSSFP_ihMT(flips(ii),b1sqrd{ii,1},Delta_Hz,TR,tau,dphi,tissuepars);
    %%% Asymmetric (2-band) case
    Mssfp(ii,2) = ssSSFP_ihMT(flips(ii),b1sqrd{ii,2},Delta_Hz,TR,tau,dphi,tissuepars);
    %%% Single band case
    Mssfp(ii,3) =ssSSFP_ihMT(flips(ii),b1sqrd{ii,1}.*[0 1 0],Delta_Hz,TR,tau,dphi,tissuepars);
end


%% SPGR steady state

%%% Extra pars needed for SSFP
dphi=0;

% loop over the different flips
Mspgr = zeros(nf,3);
Delta_Hz = [-delta 0 delta];
for ii=1:nf
    %%% Symmetric (3-band) case
    Mspgr(ii,1) = ssSPGR_ihMT(flips(ii),b1sqrd{ii,1},Delta_Hz,TR,tau,tissuepars);
    %%% Asymmetric (2-band) case
    Mspgr(ii,2) = ssSPGR_ihMT(flips(ii),b1sqrd{ii,2},Delta_Hz,TR,tau,tissuepars);
    %%% Single band case
    Mspgr(ii,3) =ssSPGR_ihMT(flips(ii),b1sqrd{ii,1}.*[0 1 0],Delta_Hz,TR,tau,tissuepars);
end

%%
figfp(21);
subplot(1,2,1)
plot(rad2deg(flips),abs(Mssfp),'linewidth',2)
hold on
plot(rad2deg(flips),abs(freeman_hill(rad2deg(flips),TR,tissuepars.free.R1,...
    tissuepars.free.R2)),'--')

legend('3-band','2-band','1-band','Freeman-Hill')
grid on
set(gca,'fontsize',12)
xlabel('flip')
ylabel('signal')
title('Predicted signals')
% ylim([0 0.08]);

subplot(1,2,2)
plot(rad2deg(flips),abs(Mspgr),'linewidth',2)
hold on
plot(rad2deg(flips),abs(ernst(rad2deg(flips),TR,tissuepars.free.R1)),'--')

legend('3-band','2-band','1-band','Ernst')
grid on
set(gca,'fontsize',12)
xlabel('flip')
ylabel('signal')
title('Predicted signals')

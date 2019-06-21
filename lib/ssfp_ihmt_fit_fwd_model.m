% function to perform forward model

function S = ssfp_ihmt_fit_fwd_model(fit_pars,sequence_pars,flag)

%%% fit_pars contains, in this order:
% [R1f R2 M0f R1sz R1D f T2s k sf]
tissuepars=struct;
tissuepars.free.R1 = fit_pars(1);
tissuepars.free.R2 = fit_pars(2);
tissuepars.semi.M0 = fit_pars(3);
tissuepars.semi.R1 = fit_pars(4);
tissuepars.semi.R1D = fit_pars(5);
tissuepars.semi.f  = fit_pars(6);
tissuepars.semi.T2 = fit_pars(7)*1e-6;%<-- fits use microseconds units
tissuepars.k  = fit_pars(8);

nf = length(sequence_pars.flips);

%%% ****** Bieri Scheffler correction ****** 
R1 = tissuepars.free.R1;
R2 = tissuepars.free.R2;

trfe = sequence_pars.tau * 1.2/sequence_pars.TBP;
zeta = 0.68 - 0.125 * (1+trfe/sequence_pars.TR)*R1/R2;
R2hat = (1 - zeta*trfe/sequence_pars.TR)*R2;

%%% Modify 
tissuepars.free.R2 = R2hat;

%%% Lineshape
if (nargin<3)||(isempty(flag))
    tissuepars.lineshape = 'SL';
else
    tissuepars.lineshape = 'Gaussian';
end


%%% Make the signals
S = zeros(nf,3); % 1-band 2-band 3-band
for ii=1:nf
    %%% Symmetric (3-band) case
    S(ii,3) = ssSSFP_ihMT(fit_pars(10)*sequence_pars.flips(ii),fit_pars(10)^2*sequence_pars.b1sqrd{ii,3},...
    sequence_pars.Delta_Hz,sequence_pars.TR,sequence_pars.tau,sequence_pars.dphi,tissuepars);
    %%% Asymmetric (2-band) case
    S(ii,2) = ssSSFP_ihMT(fit_pars(10)*sequence_pars.flips(ii),fit_pars(10)^2*sequence_pars.b1sqrd{ii,2},...
        sequence_pars.Delta_Hz,sequence_pars.TR,sequence_pars.tau,sequence_pars.dphi,tissuepars);
    %%% Single band case
    S(ii,1) =ssSSFP_ihMT(fit_pars(10)*sequence_pars.flips(ii),fit_pars(10)^2*sequence_pars.b1sqrd{ii,1},...
        sequence_pars.Delta_Hz,sequence_pars.TR,sequence_pars.tau,sequence_pars.dphi,tissuepars);
end

if fit_pars(9)==-1 %<--- this is a flag to make it normalise
    S = abs(S)/mean(abs(S(:)));
else
    % scale by 9th parameter
    S = fit_pars(9)*abs(S);
end
end
function [tissuepars,R2hat] = Bieri_Scheffler_finite_correction(TBP,tau,TR,tissuepars)
% 15-2-2019: technically only on resonance, there is an off resonance
% version
R1 = tissuepars.free.R1;
R2 = tissuepars.free.R2;

trfe = tau * 1.2/TBP;
zeta = 0.68 - 0.125 * (1+trfe/TR)*R1/R2;
R2hat = (1 - zeta*trfe/TR)*R2;

%%% Modify 
tissuepars.free.R2 = R2hat;

end

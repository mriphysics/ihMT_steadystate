%%% function [Mss,Mz] = ssSSFP_ihMT(flipangle,b1sqrd,Delta_Hz,TR,tau,dphi, tissuepars,varargin)
%
%   Steady-state ihMT SPGR sequence. For non-selective multiband sequences
%
%   INPUTS:         
%           flipangle  = flip angle on resonance (rad)
%           b1sqrd     = mean square B1+ per frequency band of
%                        multiband pulse (over the duration of the pulse). 
%                        1x3 vector (so far we assume 3 bands maximum but 
%                        this could be changed). Units uT^2
%           Delta_Hz   = 1x3 vector of frequencies of each band. Typically 
%                        [-delta 0 delta]. Units Hz
%           TR         = repetition time, sec
%           tau        = pulse duration, sec
%           dphi       = Off-resonance phase gained per TR, unit=radians
%           tissuepars = structure containing all tissue parameters. See
%                        init_tissue()
%
%  OUTPUTS:
%           Mss        = Steady-state Mxy (after excitation pulse)
%           Mz         = Longitudinal magnetization, including semisolid
%                        terms
%
% (c) Shaihan Malik 2019. King's College London

function [Mss,Mz] = ssSSFP_ihMT(flipangle,b1sqrd,Delta_Hz,TR,tau,dphi, tissuepars,varargin)

%%% Unpack tissue parameters
M0s = tissuepars.semi.M0;  %<--- TOTAL semisolid fraction
f   = tissuepars.semi.f;%<--- f = fraction of semisolid pool that is dipolar-coupled
M0f = 1-M0s;            %<--- free fraction is 1-M0s

R1  = [tissuepars.free.R1 tissuepars.semi.R1 tissuepars.semi.R1 tissuepars.semi.R1D];
R2f = tissuepars.free.R2;

% semisolid T2, used in lineshape calculation
T2s = tissuepars.semi.T2;

% Overall exchange rate for free and both semisolid pools
k = tissuepars.k;

%%% lineshape
switch tissuepars.lineshape
    case 'SL'
        [G,w_loc] = SuperLorentzian_lineshape(T2s,Delta_Hz);% seconds
    case 'Gaussian'
        [G,w_loc] = gauss_lineshape(T2s,Delta_Hz);% seconds
end

% gamma for RF calculation
gam = 267.5221; %< rad /s /uT



%% Evolution during RF pulses

Omega_semi = zeros(3);
% loop over frequency bands
for ii=1:3
    W = pi*gam^2*b1sqrd(ii)*G(ii);
    D = 2*pi*Delta_Hz(ii)/w_loc;
    Omega_semi = Omega_semi + [[-W 0 0];[0 -W W*D];[0 W*D -W*D^2]];
end


%%% Assume RF pulse has zero phase (i.e. B1x only) - the matrix exponential
%%% of this times tau is simply a rotation matrix
Omega_free = [0 0 0;0 0 flipangle/tau;0 -flipangle/tau 0];

R = blkdiag(expm(Omega_free*tau),expm(Omega_semi*tau));

    
%%% Evolution in free precesion periods

% Free pool transverse components
La = [-R2f dphi/TR;-dphi/TR -R2f];

% the rest
Lb = [-k*M0s-R1(1) k*M0f k*M0f 0;k*(1-f)*M0s -k*M0f-R1(2) 0 0;...
    k*f*M0s 0 -k*M0f-R1(3) 0;0 0 0 -R1(4)];

Lambda = blkdiag(La,Lb);

C = [0 0 R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0].';

% phase alternation matrix
Phi = diag([-1 -1 1 1 1 1]);
        
        
%% Now compute the steady-state
S = expm(Lambda*TR);
I = eye(6);

% SS solution
Mss_all =  S^0.5 * inv(Phi - R*S) * R * (S-I)*inv(Lambda)*C;

Mss = ([1 1i 0 0 0 0])*Mss_all; % return transverse magnetization

Mz = abs(Mss_all(3:end)); % return Mz



end
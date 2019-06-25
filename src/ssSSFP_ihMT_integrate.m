%%% [Mss,Mz] = ssSSFP_ihMT_integrate(b1pulse,dt,Delta_Hz,TR,dphi, tissuepars)
%
%   Steady-state ihMT bSSFP sequence with eigenvector based time integration method
%
%   INPUTS:         
%           b1pulse    = RF pulse, Mx3 array (M=#timepoints, 3=frequency
%                        bands). Units are uT
%           dt         = dwell time, sec
%           Delta_Hz   = 1x3 vector of frequencies of each band. Typically 
%                        [-delta 0 delta]. Units Hz
%           TR         = repetition time, sec
%           dphi       = off-resonance phase per TR (rad)
%           tissuepars = structure containing all tissue parameters. See
%                        init_tissue()
%
%  OUTPUTS:
%           Mss        = Steady-state Mxy (after excitation pulse)
%           Mz         = Longitudinal magnetization, including semisolid
%                        terms
%
% (c) Shaihan Malik 2019. King's College London

function [Mss,Mz] = ssSSFP_ihMT_integrate(b1pulse,dt,Delta_Hz,TR,dphi, tissuepars)

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
        [G,w_loc] = SuperLorentzian_lineshape(T2s,Delta_Hz,'interpzero');% seconds
    case 'Gaussian'
        [G,w_loc] = gauss_lineshape(T2s,Delta_Hz);% seconds
end

% gamma for RF calculation
gam = 267.5221; %< rad /s /uT



%% Lambda matrix and C are time invariant

% Free pool transverse components
La = [-R2f dphi/TR;-dphi/TR -R2f];

% the rest
Lb = [-k*M0s-R1(1) k*M0f k*M0f 0;k*(1-f)*M0s -k*M0f-R1(2) 0 0;...
    k*f*M0s 0 -k*M0f-R1(3) 0;0 0 0 -R1(4)];

Lambda = blkdiag(La,Lb);

C = [0 0 R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0].';

% Augmented matrix form - this is a rate so applies whenever
Atilde = cat(1,[Lambda C],zeros(1,7));

% Work out pulse duration, and hence the duration of the evolution period.
% We can simply apply the constant rate for that whole period with no
% approximation
nt = size(b1pulse,1); 
tau = dt*nt;

% Evolution in time between pulses (TR-tau)
Xtilde_norf = expm(Atilde*(TR-tau));

% 24-5-2019: Also re-generate the full multiband pulse
b1_multiband = zeros(nt,1);
tt = (1:nt)'*dt;
for ii=1:3
    b1_multiband = b1_multiband + exp(2*pi*Delta_Hz(ii)*1i*tt).*b1pulse(:,ii);
end

%% Now look at the RF matrices and integrate over pulse


% initialise this matrix
Xtilde_rf = eye(7);

for tt=1:nt
    
    % 24-5-2019: Apply full multiband pulse
    b1x = real(b1_multiband(tt));
    b1y = imag(b1_multiband(tt));
    OmegaFree = gam*[0 0 -b1y;0 0 b1x;b1y -b1x 0];
    
    % now loop over the bands
    OmegaSemi = zeros(3,3);

    for ii=1:3
        w1 = gam*abs(b1pulse(tt,ii));
        w = pi*w1^2*G(ii);
        D = 2*pi*Delta_Hz(ii)/w_loc;
        OmegaSemi = OmegaSemi + [[-w 0 0];[0 -w w*D];[0 w*D -w*D^2]];
    end
    
    Omega = blkdiag(OmegaFree,OmegaSemi);
    
    % Now make overall evolution matrix
    Atilde = cat(1,[(Lambda+Omega) C],zeros(1,7));
    
    % apply matrix exponential and multiply into existing matrix
    Xtilde_rf = expm(Atilde*dt)*Xtilde_rf;
    
end
    

%% Now compile these into a sequence prediction

% phase alternation matrix
Phi = diag([-1 -1 1 1 1 1 1]);

% combine all. First apply evolution, then RF pulse, then phase
% alternation. This is the steady state from the perspective of right at
% the end of each RF pulse
%X = Phi * Xtilde_rf * Xtilde_norf;

% Switch reference time to the middle of the TR
X = Phi * Xtilde_norf^0.5 *Xtilde_rf * Xtilde_norf^0.5;


%%% Take eigenvector decomposition
[v,d,w] = eig(X);
        
MssAll = v(:,end)/v(end,end); %<-- normalise to the last component since this should always be 1

% Component to return, Mx+iMy
S = [1 1i 0 0 0 0 0];

Mss = S*MssAll;
Mz = abs(MssAll(3:end-1)); % return Mz
        

end
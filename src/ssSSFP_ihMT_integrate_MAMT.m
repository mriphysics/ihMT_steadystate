%%% [Mss,Mz] = ssSSFP_ihMT_integrate(b1pulse,dt,Delta_Hz,TR,dphi,NTR, tissuepars)
%
%   Steady-state ihMT SPGR sequence using time integration with MAMT
%   (Portnoy and Stanisz MRM 2007) style approach
%
%   INPUTS:         
%           b1pulse    = RF pulse, Mx3 array (M=#timepoints, 3=frequency
%                        bands). Units are uT
%           dt         = dwell time, sec
%           Delta_Hz   = 1x3 vector of frequencies of each band. Typically 
%                        [-delta 0 delta]. Units Hz
%           TR         = repetition time, sec
%           dphi       = off-resonance phase per TR (rad)
%           NTR        = number of TR periods to simulate
%           tissuepars = structure containing all tissue parameters. See
%                        init_tissue()
%
%  OUTPUTS:
%           Mss        = Steady-state Mxy (after excitation pulse)
%           Mz         = Longitudinal magnetization, including semisolid
%                        terms
%           Mt         = time resolved calculation from end of each TR
%
% (c) Shaihan Malik 2019. King's College London

function [Mss,Mz,Mt] = ssSSFP_ihMT_integrate_MAMT(b1pulse,dt,Delta_Hz,TR,dphi,NTR, tissuepars)

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

if nargout==3 %<--- 3rd output is time trace
    timeresolved=true;
    Mt = zeros(6,NTR);
else
    timeresolved=false;
end


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

% Start at equilibrium
M = [0 0 M0f (1-f)*M0s f*M0s 0].';

% initialise identity matrix
I = eye(6);

% phase alternation matrix
Phi = diag([-1 -1 1 1 1 1]);

% These can be defined outside the loop
eL = expm(Lambda*(TR-tau)); %<--- could be defined outside loop
eL_con = (eL-I)*inv(Lambda)*C;
        
% Loop over TR periods
for ix = 1:NTR

    % Right at the start, apply phase alternation matrix
    M = Phi*M;
    
    % Now loop over samples of RF Pulse
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
        A = Lambda+Omega;
        
        % apply forward solution as per MAMT paper
        eA = expm(A*dt);
        M = eA*M + (eA-I)*inv(A)*C;
        
    end %<-- end of RF pulse

    if ix<NTR
        % Now apply the rest of the TR - only Lambda applies here, not Omega
        M = eL*M + eL_con;
    else
        % assume ix==NTR here, so it's the last run so just have half the
        % TR period
        eL = eL^0.5;
        M = eL*M + (eL-I)*inv(Lambda)*C;
    end
    if timeresolved
        Mt(:,ix)=M;
    end
end
    
        
% Component to return, Mx+iMy
S = [1 1i 0 0 0 0];

Mss = S*M;
Mz = abs(M(3:end)); % return Mz
        

end
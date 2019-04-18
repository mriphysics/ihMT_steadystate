%% 13-6-17: steady state with MT, balanced -- now matrix version, not analytic
%% 12-1-18: ihMT version - see lab book #22, this date
%%% Multi-band pulse with bands at delta, B1rms^2 in each band is in b1sqrd
% alpha, b1sqrd, delta all need to be nd length vectors
% 14-2-18: All B1 units are uT, all times are seconds, all rates are s^-1
% 22-10-18: this version has 2 semisolid pools - see labbook/evernote for
% this date
%   *** New parameter definitions
%       K = overall exchange constant. K = kfs/M0s=ksf/M0f
%       M0 = vector of relative M0 components = [M0f, M0s1, M0s2] - where
%       M0s1 is the NON-dipolar coupled Zeeman pool and M0s2 is the dipolar 
%       coupled Zeeman pool
%       T1 = T1 times [T1free T1semi1 T1semi2 T1dipolar] - if specified as
%       a 3-vector, assume we need to replicate the semisolid Zeeman T1
%       pools have the same Zeeman T1 
%
%       12-11-2018: added 'augmented' matrix + eigenvector approach to
%       steady state - use tag 'eig'
%
%       15-11-2018: - this version does pulse integration for computing the RF
%       transition matrix. 
%       b1pulse is an nt x 3 matrix of the pulse B1 at each band
%       this version ALWAYS uses eigenval approach
%
%       21-11-2018: - an SPGR version modified from the SSFP code. Note
%       that this IS dependent on T2 because you get relaxation during
%       pulses
function [Mss,Mz] = ssSPGR_ihMT_2pool_integrate(alpha,b1pulse, TR,trf,dphi, T1,T2,T2b,M0,K,delta,varargin)

M0f = M0(1);
M0s = M0(2)+M0(3);  %<--- TOTAL semisolid fraction
f = M0(3)/M0s;      %<--- f = fraction of semisolid pool that is dipolar-coupled

if numel(T1)==3
    T1 = [T1(1:2) T1(2) T1(3)]; % replicate the Semisolid Zeeman T1 
end
R1 = 1./T1;
R2 = 1/T2;

%%% Multi-band pulse with bands at delta, B1rms^2 in each band is in b1sqrd
nd = size(b1pulse,2);% b1sqrd, delta all need to be nd length vectors

%%% Other args
lineshape = 1; %<-- Super Lorentzian is assumed

if ~isempty(varargin)
    if strcmpi(varargin{1},'gaussian')
        lineshape = 2;
    end
      
end

%%% lineshape
switch lineshape
    case 1
        [G,w_loc] = SuperLorentzian_lineshape(T2b,delta);% seconds
    case 2
        [G,w_loc] = gauss_lineshape(T2b,delta);% seconds
end

% gamma for RF calculation
gam = 267.5221; %< rad /s /uT


%%% **** First, relaxation and exchange ****

E = diag([-R2 -R2 -R1]);
KK = blkdiag(0,0,[[-K*M0s K*M0f K*M0f];[K*(1-f)*M0s -K*M0f 0];[K*f*M0s 0 -K*M0f]],0);

% for SPGR can ignore dephasing - we later assume perfect spoiling
%gammaB0  = diag([-1i*dphi/TR 1i*dphi/TR 0 0 0 0]); % phi is defined as dephasing per TR

% Lambda = E+KK+gammaB0;
Lambda = E+KK;
Z = [0 0 R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0].';

% Augmented matrix form - this is a rate so applies whenever
LambdaPrime = cat(1,[Lambda Z],zeros(1,7));


% Evolution in time between pulses (TR-trf)
expLambdaPrime = expm(LambdaPrime*(TR-trf));


%%% *** Now look at the RF matrices and integrate over pulse


dt = 6.4e-6; % should input this

nt = size(b1pulse,1); %

% initialise this matrix
expRF = eye(7);

for tt=1:nt
    
    % something like this:
    w1 = gam*b1pulse(tt,2);% central band
    OmegaFree = [[0 0 1i*w1];[0 0 -1i*conj(w1)];[1i*conj(w1)/2 -1i*w1/2 0]];
    
    % now loop over the bands
    OmegaSemi = zeros(3,3);

    for ii=1:nd
        w1 = gam*abs(b1pulse(tt,ii));
        w = pi*w1^2*G(ii);
        D = 2*pi*delta(ii)/w_loc;
        OmegaSemi = OmegaSemi + [[-w 0 0];[0 -w w*D];[0 w*D -w*D^2]];
    end
    
    Omega = blkdiag(OmegaFree,OmegaSemi);
    
    % Now make overall evolution matrix
    A = cat(1,[(Lambda+Omega) Z],zeros(1,7));
    
    % apply matrix exponential and multiply into existing matrix
    expRF = expm(A*dt)*expRF;
    
end
    

%% Now compile these into a sequence prediction
% 21-11-2018: this part is different for SPGR, above is all the same

% Spoiling matrix
D = diag([0 0 1 1 1 1 1]); %<-- kill transverse components during evolution period

% combine all
X = expRF * D * expLambdaPrime; %<--- this is directly after pulse, since expRF is applied LAST (left-most)

%%% Take eigenvector decomposition
[v,d,w] = eig(X);
        
MssAll = v(:,end)/v(end,end); %<-- normalise to the last component so it stays as 1

% Component to return = M+
S = [1 0 0 0 0 0 0];

Mss = S*MssAll;
Mz = abs(MssAll(3:end-1)); % return Mz
        

end
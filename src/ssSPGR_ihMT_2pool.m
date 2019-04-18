%% Steady state SPGR function for dipolar relaxation 21/12/2017]
% initial notes in lab book 21/12/2017
% 7-2-2018: assume all flips are the same
% 7-2-2018: rewrite using alternative notation and only consider dipolar T1
% during pulse - see lab book
% 13-2-18: All B1 units are uT, all times are seconds, all rates are s^-1
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
function [Mss,Mss_all] = ssSPGR_ihMT_2pool(alpha,b1sqrd,T1,T2b,M0,K,delta,TR,trf,varargin)

M0f = M0(1);
M0s = M0(2)+M0(3);  %<--- TOTAL semisolid fraction
f = M0(3)/M0s;      %<--- f = fraction of semisolid pool that is dipolar-coupled

if numel(T1)==3
    T1 = [T1(1:2) T1(2) T1(3)]; % replicate the Semisolid Zeeman T1 
end
R1 = 1./T1; %<--- inverse seconds

%%% Convert k values to inverse seconds (input is inverse ms)


%%% Multi-band pulse with bands at delta, B1rms^2 in each band is in b1sqrd
nd = length(b1sqrd);% b1sqrd, delta all need to be nd length vectors

%%% Other args
lineshape = 1; %<-- Super Lorentzian is assumed
version = 1; %<-- standard SS formalism vs eigenvector version
if ~isempty(varargin)
    if strcmpi(varargin{1},'gaussian')
        lineshape = 2;
    end
    
    if strcmpi(varargin{1},'eig')
        version = 2;
    end
    
end


%%% lineshape
switch lineshape
    case 1
        [G,w_loc] = SuperLorentzian_lineshape(T2b,delta);% seconds
    case 2
        [G,w_loc] = gauss_lineshape(T2b,delta);% seconds
end

%%% Build up RF interaction matrix for Zeeman & Dipolar interaction
gam = 267.5221; %< rad /s /uT
W = zeros(3);
% loop over frequency bands
for ii=1:nd
    w = pi*gam^2*b1sqrd(ii)*G(ii); 
    D = 2*pi*delta(ii)/w_loc;
    W = W + [[-w 0 0];[0 -w w*D];[0 w*D -w*D^2]];
end

Alpha = blkdiag(cos(alpha),expm(W*trf));


%%% Evolution in free precesion periods
E = diag(-R1);
KK = blkdiag([[-K*M0s K*M0f K*M0f];[K*(1-f)*M0s -K*M0f 0];[K*f*M0s 0 -K*M0f]],0);
Lambda = E+KK;
Z = [R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0].';

% Component to return
S = [sin(alpha) 0 0 0];
        
%% Now decide on which version:

switch version
    
    case 1  % Steady-state formula version
        Xi = expm(Lambda*(TR));
        I = eye(4);
        Phi = (Xi-I)*inv(Lambda)*Z;
                       
        % SS solution
        Mss_all =  inv(I - Alpha*Xi) * (Alpha * Phi);
        Mss = S*Mss_all; % return signal

    case 2  % Eigendecomposition version

        LC = cat(1,[Lambda Z],zeros(1,5));
        
        %%% RF matrix
        AlphaPrime = blkdiag(Alpha,1); %<-- add 1 here not zero as this is post expm
        
        %%% overall matrix
        X = expm(LC*TR) * AlphaPrime; %<-- want state BEFORE RF pulse, as we multiple by sin theta
        
        [v,d,w] = eig(X);
        
        MssAll = v(1:4,5)/v(5,5); %<-- normalise to the last component so it stays as 1 .... dubious?

        Mss = S*MssAll;

end
end
%%% New function, Feb 2019 for publication


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

%%% Other args
% lineshape = 1;  %<-- Super Lorentzian is assumed
version = 1;    %<-- standard SS formalism vs eigenvector version

if ~isempty(varargin)
%     if strcmpi(varargin{1},'gaussian')
%         lineshape = 2;
%     end
    
    if strcmpi(varargin{1},'eig')
        version = 2;
    end
    
end

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
        
        
%% Now decide on which version:
switch version
    
    case 1  % Steady-state formula version
        S = expm(Lambda*TR);
        I = eye(6);
            
        % SS solution
        Mss_all =  S^0.5 * inv(Phi - R*S) * R * (S-I)*inv(Lambda)*C;
                
        Mss = ([1 1i 0 0 0 0])*Mss_all; % return transverse magnetization
        
        Mz = abs(Mss_all(3:end)); % return Mz

    case 2 % Eigenvector decomposition version
%         LC = cat(1,[Lambda Z],zeros(1,7));
%         
%         %%% RF matrix
%         AlphaPrime = blkdiag(Alpha,1); %<-- add 1 here not zero as this is post expm
%         
%         %%% spoiling
%         D = blkdiag(D,1);
%         
%         %%% overall matrix
%         X = AlphaPrime * D * expm(LC*TR);
%         
%         [v,d,w] = eig(X);
%         
%         MssAll = v(1:end-1,end)/v(end,end); %<-- normalise to the last component so it stays as 1 .... dubious?
% 
%         Mss = S*MssAll;
%         Mz = abs(MssAll(3:end)); % return Mz
end

end
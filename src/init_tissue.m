%%% Feb 2019 Initialise tissue parameters
function pars = init_tissue(name)

pars = struct;

if nargin==0
    % set to default
    name = 'ic';
end
    
% Replace all letters with lower case
name = lower(name);

% initialise the empty structure
% Rates are all s^-1, times on seconds, magnetizations referenced to total
% M0=1

%%% Free water pool
pars.free.R1 = [];
pars.free.R2 = [];
%%% semisolid pool
pars.semi.M0 = [];   % assume M0f+M0s = 1
pars.semi.R1 = [];   %<--- assume same for both s1 and s2
pars.semi.R1D = [];
pars.semi.f  = [];   % fraction of M0s that is in s1
pars.semi.T2 = [];   %<--- assume same for both s1 and s2
%%% overall exchange rate constant
pars.k  = [];

switch name
    
    case 'simple'
                
        pars.free.R1 = 1;
        pars.free.R2 = 10; % 100ms T2
        
        pars.semi.M0 = 0.2;    
        pars.semi.R1 = 1;      
        pars.semi.R1D = 1/10e-3;
        pars.semi.f  = 0.65;    
        pars.semi.T2 = 12e-6;   
        
        pars.k  = 50;

        pars.lineshape = 'SL'; %<-- super Lorentzian
        
    case 'ic'
        
        %%% Parameters from Mchinda 2017 -- Default parameters
        
        pars.free.R1 = 1/650e-3;
        pars.free.R2 = 1/80e-3;
        
        %pars.semi.M0 = 0.17;    %<--- assume that total M0 is 1
        pars.semi.M0 = 0.147;    %<--- M0f=1 for Varma work. Here M0s = M0s_Varma / (1+M0s_Varma)
        pars.semi.R1 = 1;       %<--- assume same for both s1 and s2
        pars.semi.R1D = 1/6.5e-3; %<--- dipolar relaxation rate
        pars.semi.f  = 0.65;    %<--- fraction that has long T1D
        pars.semi.T2 = 12e-6;   %<--- assume same for both s1 and s2
        
        %%% overall exchange constant
        pars.k  = 65;

        pars.lineshape = 'SL';

end

end
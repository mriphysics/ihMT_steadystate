% Function to compute Super Lorentzian lineshape and associated local field
% strength, w_loc.
% 
% [G,w_loc] = SuperLorentzian_lineshape(T2b,fsample,varargin)
% 
% INPUTS:       T2b = T2 of semisolid compartment, seconds
%               fsample   = frequencis at which function is to be evaluated (can
%               be vector)
%               optional: 'interpzero' - interpolate between ±1.5kHz as per
%               Gloor et al 2008, to remove high peak at f=0.
%
% OUTPUTS:      G   = lineshape value, s
%               w_loc = local field term (rad/s)
%               
%               optional:
% Shaihan Malik (c), King's College London, April 2019

function [G,w_loc] = SuperLorentzian_lineshape(T2b,fsample,varargin)

%%% A suggested by Gloor, we can interpolate the lineshape across from
%%% ±1kHz
interpzero = false;
for ii=1:length(varargin)
    if strcmpi(varargin{ii},'interpzero')
        interpzero = true;
    end
end


if interpzero
    % compute over a wider range of frequencies
    n=128;
    ff = linspace(1.1*(min(fsample)),1.1*(max(fsample)),n);
else
    ff = fsample;
end

%%% Variables for SL, predefine
th = linspace(0,pi/2,512);
dth = th(2)-th(1);

[thg,ffg]=meshgrid(th,ff);

g = sin(thg).*sqrt(2/pi).*T2b./(abs(3*cos(thg).^2-1));
g = g .* exp(-2*(2*pi*ffg*T2b./abs(3*cos(thg).^2-1)).^2);
G = dth*sum(g,2);

%
if interpzero
    po = find(abs(ff)<1.5e3); % points to interpolate
    pu = find((abs(ff)>1e3)&(abs(ff)<2e3)); % points to use
    
    Gi = spline(ff(pu),G(pu),ff(po));
    G(po) = Gi;
    
    %%% now finally interpolate the required frequencies
    G = interp1(ff,G,fsample);
end

w_loc = 1/(sqrt(15)*T2b);


end



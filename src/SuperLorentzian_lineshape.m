% Function to compute Super Lorentzian lineshape and associated local field
% strength, w_loc.
% 
% [G,w_loc] = SuperLorentzian_lineshape(T2s,fsample,varargin)
% 
% INPUTS:       T2s = T2 of semisolid compartment, seconds
%               fsample   = frequencis at which function is to be evaluated (can
%               be vector)
%               optional: 'interpzero' - interpolate between ±1.5kHz as per
%               Gloor et al 2008, to remove high peak at f=0.
%
% OUTPUTS:      G   = lineshape value, s
%               w_loc = local field term (rad/s)
%               
%               
% Shaihan Malik (c), King's College London, April 2019

function [G,w_loc] = SuperLorentzian_lineshape(T2s,fsample,varargin)

%%% A suggested by Gloor, we can interpolate the lineshape across from
%%% ±1kHz
interpzero = false;
interptheta = false;
Ntheta=128; %<-- number of points for theta integral
for ii=1:length(varargin)
    if strcmpi(varargin{ii},'interpzero')
        interpzero = true;
    end
    if strcmpi(varargin{ii},'Ntheta')
        Ntheta = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'interptheta')
        interptheta = true;
    end
end


if interpzero
    % compute over a wider range of frequencies
    n=128;
    if min(fsample)>-2e3
        fmin=-2e3;
    else
        fmin = min(fsample)*1.1;
    end
    if max(fsample)<2e3
        fmax=2e3;
    else
        fmax = max(fsample)*1.1;
    end
    ff = linspace(fmin,fmax,n);
else
    ff = fsample;
end

%%% Variables for SL, predefine
th = linspace(0,pi/2,Ntheta);
dth = th(2)-th(1);

[thg,ffg]=meshgrid(th,ff);

g = sin(thg).*sqrt(2/pi).*T2s./(abs(3*cos(thg).^2-1));
g = g .* exp(-2*(2*pi*ffg*T2s./abs(3*cos(thg).^2-1)).^2);

%% at this point try to interpolate over the angles not the line value
if interptheta
   
    magic_angle = acos(1/sqrt(3));

    po = find(abs(th-magic_angle)<deg2rad(3)); % points to interpolate
    pu = find((abs(th-magic_angle)>deg2rad(3))&(abs(th-magic_angle)<deg2rad(10))); % points to use
   
   %%% loop over frequencies
   for ii=1:length(ff)
       tmp = spline(th(pu),g(ii,pu),th(po));
       g(ii,po) = tmp;
   end
end

%%
G = dth*sum(g,2);

%
if interpzero
    po = find(abs(ff)<1.e3); % points to interpolate
    pu = find((abs(ff)>1.e3)&(abs(ff)<2e3)); % points to use
    
    Gi = spline(ff(pu),G(pu),ff(po));
    G(po) = Gi;
    
    %%% now finally interpolate the required frequencies
    G = interp1(ff,G,fsample);
end

w_loc = 1/(sqrt(15)*T2s);


end



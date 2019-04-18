%% Super Lorentzian lineshape including second moment, and also interpolated.
% return LUT interpolated from over range of the specified sample
% frequencies
% Units of G are seconds 
% 7-2-2018

function [G,D,ff] = SuperLorentzian_lineshape(T2b,fsample,varargin)

%%% A suggested by Gloor, we can interpolate the lineshape across from
%%% ±1kHz
interpzero = false;
for ii=1:length(varargin)
    if strcmpi(varargin{ii},'interpzero')
        interpzero = true;
    end
end
interpzero = true;
%%% define frequency range
% n=512;
n=128;
%ff = linspace(-30e3,30e3,n);
ff = linspace(1.1*(min(fsample)),1.1*(max(fsample)),n);

%%% compute G for this range
G = zeros([n 1]);

%%% Variables for SL, predefine
th = linspace(0,pi/2,128);
dth = th(2)-th(1);

% tic
% for ii=1:n
%     G(ii) = SL(ff(ii));
% end
% toc
%% Try to do as matrix

[thg,ffg]=meshgrid(th,ff);

g = sin(thg).*sqrt(2/pi).*T2b./(abs(3*cos(thg).^2-1));
g = g .* exp(-2*(2*pi*ffg*T2b./abs(3*cos(thg).^2-1)).^2);
G = dth*sum(g,2);

%%
if interpzero
    %disp('SL:interp')
    %%% interpolate
    po = find(abs(ff)<1.5e3); % points to interpolate
    pu = find((abs(ff)>1e3)&(abs(ff)<2e3)); % points to use
    
    Gi = spline(ff(pu),G(pu),ff(po));
    G(po) = Gi;
end

D = 1/(sqrt(15)*T2b);


%%% If the user requested samples, now need to sample the correct values
if nargin>1
    G = interp1(ff,G,fsample);
end
    
%     function gg = SL(f)
%         g = sin(th).*sqrt(2/pi).*T2b./(abs(3*cos(th).^2-1));
%         g = g .* exp(-2*(2*pi*f*T2b./abs(3*cos(th).^2-1)).^2);
%         gg = dth*sum(g);
%     end
end



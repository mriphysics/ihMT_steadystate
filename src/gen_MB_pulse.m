function [pulseMB, b1sqrd, pulse_per_band,TBP] = gen_MB_pulse(theta,tau,TR,b1rms_total,delta,nband,varargin)
%%% function based on publication
%%% [pulseMB, b1sqrd, pulse_per_band,TBP] = gen_MB_pulse(theta,tau,TR,b1rms_total,delta,nband,varargin)
%%%   	arguments:
% %   	-  theta: 		on-resonance flip angle (rad)
% % 	-  tau:   		pulse duration (s)
% % 	-  TR:    		repetition time (s)
% % 	-  b1rms_total: overall RMS B1 of the whole sequence (uT)
% % 	-  delta:	    offset frequency for off-resonant bands (Hz) 
% % 	-  nband:		string argument for number of bands. Options are '2+', '2-', or '3'.
% % 	
% % 	optional arguments:
% % 	-  alpha:		scalar value determining width of Gaussian envelope for pulse (default is 3)
% % 	-  dt:			RF sampling duration. Default 6.4us
%%% Feb 2019
%%% April 2009 - output the TBP as well
%%% 
%%
gausswin_alpha = 3; %<-- shape of basic pulse
dt = 6.4e-6;

if ~isempty(varargin)
    for ii=1:length(varargin)
        if strcmpi(varargin{ii},'alpha')
            gausswin_alpha = varargin{ii+1};
        end
        
        if strcmpi(varargin{ii},'dt')
            dt = varargin{ii+1};
        end
        
        if ~isstr(nband)
            disp('gen_MB_pulse: warning, nband should be a string, either 2+ 2- or 3 - converting...')
            if nband==2
                disp('gen_MB_pulse: assuming 2+')
                nband = '2+';
            else
                if nband==3
                    nband = '3';
                end
            end
        end
    end
end

%% set up time domain
nt = ceil(tau/dt);
tau = dt*nt;

%% basic shape
pulse = gausswin(nt,gausswin_alpha);
pulse = pulse - min(pulse);
pulse = pulse / max(pulse);


gam = 267.5221; %< rad /s /uT
teff =  (gam*sum(pulse)*dt)/(gam*max(pulse))/tau; %<-- teff of shape - normalised to duratio


%% Make pulse 

tt = dt*(1:nt);

b1_max_ONres = (theta)/(gam*teff*tau); % uT units

%%% Create on resonance pulse
pulseSB = pulse * b1_max_ONres;

%%% compute RMS B1 of on resonance pulse
b1_rms_ONres = sqrt(sum(abs(pulseSB).^2)*dt/TR);

%%% If this is already higher than specified total B1rms, then a multiband
%%% solution is not possible
if (b1_rms_ONres > b1rms_total)
    disp('WARNING: Single band already exceeds specified total B1RMS - no solution is possible');
    return
end

%%% compute beta factor - ratio of off resonant to on resonant power
beta = (b1rms_total^2-b1_rms_ONres^2)/b1_rms_ONres^2;

%%% Modulation functions
switch nband
    case '2+'
        wt = 1 + sqrt(beta)*(cos(2*pi*delta*tt) + 1i*sin(2*pi*delta*tt));
        
    case '2-'
        wt = 1 + sqrt(beta)*(cos(2*pi*delta*tt) - 1i*sin(2*pi*delta*tt));
        
    case '3'
        wt = 1 + sqrt(2*beta)*cos(2*pi*delta*tt);
        
    otherwise
        disp('gen_MB_pulse: error, number of bands incorrectly defined')
        return
end


%%% Make modulated pulse
pulseMB = pulseSB(:) .* wt(:);

%%% Maximum
b1_max_pulse = max(abs(pulseMB));


fprintf(1,'Flip = %1.1f, B1 max = %1.3f uT, B1 rms = %1.3f uT\n',r2d(gam*dt*abs(sum(pulseMB))),b1_max_pulse,sqrt(sum(abs(pulseMB).^2)*dt/TR))

%%% Work out power in each band - not too hard
switch nband
    case '2+'
        b1sqrd = b1rms_total^2 * TR/tau * [0 1 beta]/(1+beta);
        
    case '2-'
        b1sqrd = b1rms_total^2 * TR/tau * [beta 1 0]/(1+beta);
        
    case '3'
        b1sqrd = b1rms_total^2 * TR/tau * [beta/2 1 beta/2]/(1+beta);
        
end



if nargout>2
    % we need to also put out the individual pulses per band - don't need
    % phase shift here as it is only used for B1^2
    pulse_per_band = zeros([nt 3]);
    pulse_per_band(:,2) = pulseSB;   
    
    switch nband
        case '2+'
            pulse_per_band(:,3) = pulseSB*sqrt(beta);
            
        case '2-'
            pulse_per_band(:,1) = pulseSB*sqrt(beta);
            
        case '3'
            pulse_per_band(:,1) = pulseSB*sqrt(beta/2);
            pulse_per_band(:,3) = pulseSB*sqrt(beta/2);
    end

end

%%% Finally TBP - based on Eq.A2 from Bieri & Scheffler paper
stdev = (nt-1)/(2*gausswin_alpha);%<--- from documentation of 'gausswin'
sig = stdev/(nt) * tau;% Scale to time units
fwhm = sqrt(2*log(2))/(pi*sig);%<--- Eq.A2
TBP = fwhm * tau;

end
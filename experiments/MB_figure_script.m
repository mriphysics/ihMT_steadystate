%% 19-4-2019: This script makes Figure 3 from the paper
% Shaihan Malik, King's College London, 2019


theta = deg2rad(35);
dur = 2e-3;
TR = 5e-3;
b1rms = 2;
delta = 7e3;
nband = 2;
dt=6e-6;
gauss_alpha = 3;

pulses={};
pulses{1} = gen_MB_pulse(theta,dur,TR,b1rms,delta,'2+','alpha',gauss_alpha,'dt',dt);
pulses{2} = gen_MB_pulse(theta,dur,TR,b1rms,delta,'2-','alpha',gauss_alpha,'dt',dt);
pulses{3} = gen_MB_pulse(theta,dur,TR,b1rms,delta,'3','alpha',gauss_alpha,'dt',dt);



%% Examine RF power in spectral bands
ff = linspace(-10e3,10e3,1000)';
df=ff(2)-ff(1);
nt = length(pulses{1});%<- all same length
tt = dt*(1:nt);
F = exp(-1i*2*pi*ff*tt)*(dt*1e3)/sqrt(numel(ff)); %<- time in ms

b1sqrd_tau = {}; % integral of B1^2 dt, units uT^2 ms
b1sqrd = {}; % average of above, over pulse duration

% frequency bands to look in
band_ix = {};
bw=2e3;
band_ix{1} = find((ff>-(delta+bw/2))&(ff<-(delta-bw/2)));
band_ix{2} = find((ff>-bw/2)&(ff<bw/2));
band_ix{3} = find((ff>(delta-bw/2))&(ff<(delta+bw/2)));

for ii=1:3
    
    % power spectrum in F
    pwr_spec = abs(F*pulses{ii}).^2;
    
    % look at power in each band - this is scaled to have units
    % equivalent to B1^2 in uT^2
    b1sqrd_tau{ii} = zeros([1 3]);
    for kk=1:3
        b1sqrd_tau{ii}(kk) = sum(pwr_spec(band_ix{kk}))*df;
    end
    tau = 1e3*dt*nt; %<-- this is constant
    b1sqrd{ii} = b1sqrd_tau{ii}/tau;
    
end
        


%% Figure with time domain as well

figfp(202);

tag = {'2^+','2^-','3'};
for ii=1:3
    pwr = abs(F*pulses{ii}).^2;
    
    subplot(3,2,ii*2-1)
    plot(ff*1e-3,pwr,'linewidth',1,'color',[0 0 0])
    grid on
    ylim([0 0.025])
    xlabel('\Delta/2\pi (Hz)')
    ylabel('power (au)')
    set(gca,'YTickLabel','','XTick',-10:2:10)
    title(sprintf('%s bands, Power spectrum',tag{ii}))
    
    %%% shadse bands
    idx = ff>-1e3&ff<1e3;
    hold on
    H = area(ff(idx)*1e-3,pwr(idx));
    set(H(1),'FaceColor',[0.7 1 0.7]);
    
    idx = ff>6e3&ff<8e3;
    hold on
    H = area(ff(idx)*1e-3,pwr(idx));
    set(H(1),'FaceColor',[1 0.5 0.5]);
    
    idx = ff>-8e3&ff<-6e3;
    hold on
    H = area(ff(idx)*1e-3,pwr(idx));
    set(H(1),'FaceColor',[1 0.5 0.5]);
    
    %%% Plot pulses
    subplot(3,2,ii*2)
    plot(tt*1e3,real(pulses{ii}))
    hold on
    plot(tt*1e3,imag(pulses{ii}))
    grid on
    xlabel('\tau, ms')
    ylabel('B_1, \muT')
    xlim([0 2])
    title('time domain')
    set(gca,'XTick',0:0.5:2,'YTick',-6:2:10)
    legend('real','imag','location','southeast')
end

%%% add annotations
gg=get(gcf,'children');
axes(gg(3));
gg(3).Position = [0.07 0.10 0.5 0.21];
text(-0.2,0.002,'1')
text(6.7,0.0025,'$\frac{\beta}{2}$','interpreter','latex','fontsize',13)
text(-7.3,0.0025,'$\frac{\beta}{2}$','interpreter','latex','fontsize',13)

axes(gg(6));
gg(6).Position = gg(3).Position + [0 0.32 0 0];
text(-0.2,0.002,'1')
text(-7.3,0.004,'$\beta$','interpreter','latex','fontsize',13)

axes(gg(9));
gg(9).Position = gg(3).Position + [0 0.64 0 0];
text(-0.2,0.002,'1')
text(6.7,0.004,'$\beta$','interpreter','latex','fontsize',13)

gg(2).Position = [0.68 0.1 0.25 0.21];
gg(5).Position = gg(2).Position + [0 0.32 0 0];
gg(8).Position = gg(2).Position + [0 0.64 0 0];
setpospap([400 200 650 600])
% print('-dpng','-r300','figs/mb_fig.png')

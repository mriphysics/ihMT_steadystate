%% Simple example for bSSFP and SPGR

tissuepars = init_tissue('simple');


%%% Set up sequence

%%% generate some pulses
flips = d2r([1:2:30 35:5:80]);
nf = length(flips);
tau = 2e-3;
TR = 5e-3;
b1_rms = 5;
delta = 7e3;
Delta_Hz = [-delta 0 delta];
dt = 10e-6;

pulses = {};
b1sqrd = {};
b1pulse= {};
for ii=1:nf
    % 3 bands
    [pulses{ii,1},b1sqrd{ii,1},b1pulse{ii,1}] = gen_MB_pulse(flips(ii),tau,TR,b1_rms,delta,'3','sigma',2,'dt',dt); 
    % 2 bands
    [pulses{ii,2},b1sqrd{ii,2},b1pulse{ii,2}] = gen_MB_pulse(flips(ii),tau,TR,b1_rms,delta,'2+','sigma',2,'dt',dt);
end


%% Simulate for both  sequences 

f = [0 1];

Mssfp = zeros(nf,3,2); % FA, #bands, f values
Mspgr = zeros(nf,3,2);

for IX = 1:2
    
    tmp = tissuepars;
    tmp.semi.f = f(IX);
    
    %%% ====== bSSFP ==========
    dphi=0;
    for ii=1:nf
        %%% Symmetric (3-band) case
        Mssfp(ii,3,IX) = ssSSFP_ihMT(flips(ii),b1sqrd{ii,1},Delta_Hz,TR,tau,dphi,tmp);
        %%% Asymmetric (2-band) case
        Mssfp(ii,2,IX) = ssSSFP_ihMT(flips(ii),b1sqrd{ii,2},Delta_Hz,TR,tau,dphi,tmp);
        %%% Single band case
        Mssfp(ii,1,IX) =ssSSFP_ihMT(flips(ii),b1sqrd{ii,1}.*[0 1 0],Delta_Hz,TR,tau,dphi,tmp);
    end
    
    %%% ====== SPGR ==========
    for ii=1:nf
        %%% Symmetric (3-band) case
        Mspgr(ii,3,IX) = ssSPGR_ihMT(flips(ii),b1sqrd{ii,1},Delta_Hz,TR,tau,tmp);
        %%% Asymmetric (2-band) case
        Mspgr(ii,2,IX) = ssSPGR_ihMT(flips(ii),b1sqrd{ii,2},Delta_Hz,TR,tau,tmp);
        %%% Single band case
        Mspgr(ii,1,IX) =ssSPGR_ihMT(flips(ii),b1sqrd{ii,1}.*[0 1 0],Delta_Hz,TR,tau,tmp);
    end
    
end

%% Plots

figfp(1)
nr=2;nc=2;

for jj=1:2
    
    subplot(nr,nc,1+(jj-1)*2)
    plot(r2d(flips),abs(Mssfp(:,:,jj)),'linewidth',1.5);
    grid on
    xlabel('Flip angle, deg')
    ylabel('Signal / M_0')
    pl = get(gca,'Children');
    pl(3).Color = [0 0 0];
    pl(1).LineStyle = '-.';
    pl(1).Color = [0.2 0.4 1];
    hold on

    
    subplot(nr,nc,2+(jj-1)*2)
    plot(r2d(flips),abs(Mspgr(:,:,jj)),'linewidth',1.5)
    grid on
    xlabel('Flip angle, deg')
    ylabel('Signal / M_0')
    
    pl = get(gca,'Children');
    pl(3).Color = [0 0 0];
    pl(1).LineStyle = '-.';
    pl(1).Color = [0.2 0.4 1];
    hold on
    
end

% add legend;
ll = legend('Single band RF pulses','2+ band RF pulses','3 band RF pulses');
ll.AutoUpdate = 'off';

%%% Now fill in some areas
gg = get(gcf,'children');

%%% For first plot fill in classic MT difference
axes(gg(5))
inBetween = [abs(Mssfp(:,1,1)); flip(abs(Mssfp(:,2,1)),1)];
x2 = [r2d(flips) flip(r2d(flips),2)];
fh=fill(x2, inBetween, 'b');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.1;

axes(gg(4))
inBetween = [abs(Mspgr(:,1,1)); flip(abs(Mspgr(:,2,1)),1)];
x2 = [r2d(flips) flip(r2d(flips),2)];
fh=fill(x2, inBetween, 'b');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.1;

axes(gg(3))
inBetween = [abs(Mssfp(:,2,2)); flip(abs(Mssfp(:,3,2)),1)];
x2 = [r2d(flips) flip(r2d(flips),2)];
fh=fill(x2, inBetween, 'r');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.1;

axes(gg(2))
inBetween = [abs(Mspgr(:,2,2)); flip(abs(Mspgr(:,3,2)),1)];
x2 = [r2d(flips) flip(r2d(flips),2)];
fh=fill(x2, inBetween, 'r');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.1;

%%% text labels
axes(gg(5))
text(30,0.11,'bSSFP','fontsize',16,'fontweight','bold')
text(-25,0.04,'f = 0','fontsize',18,'fontweight','bold','rotation',90,'fontangle','italic')

axes(gg(4))
text(30,0.055,'SPGR','fontsize',16,'fontweight','bold')

axes(gg(3))
text(-25,0.04,'f = 1','fontsize',18,'fontweight','bold','rotation',90,'fontangle','italic')

%%% annotation
a1 = annotation('textarrow',[0.2786 0.2214],[0.6878 0.7715],'String','\DeltaMT','fontsize',13);
a2 = annotation('textarrow',[0.3285 0.2713],[0.2104 0.2941],'String','\DeltaihMT','fontsize',13);
a1.HeadStyle = 'cback1';
a2.HeadStyle = 'cback1';
a1.HeadWidth = 6;
a2.HeadWidth = 6;
%
print -dpng -r300 figs/simplefigure.png
%% 19-4-2019: This script performs simulations of bSSFP and SPGR for white matter
% Shaihan Malik, King's College London, 2019

%% Generates figures 3 and 7 from the paper

tissuepars = init_tissue('ic');


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
gauss_alpha = 3;

pulses = {};
b1sqrd = {};
b1pulse= {};
for ii=1:nf
    % 3 bands
    [pulses{ii,1},b1sqrd{ii,1},b1pulse{ii,1},TBP] = gen_MB_pulse(flips(ii),tau,TR,b1_rms,delta,'3','alpha',gauss_alpha,'dt',dt); 
    % 2 bands
    [pulses{ii,2},b1sqrd{ii,2},b1pulse{ii,2},TBP] = gen_MB_pulse(flips(ii),tau,TR,b1_rms,delta,'2+','alpha',gauss_alpha,'dt',dt);
end


%% Simulate for both  sequences 

f = [0 1];

Mssfp = zeros(nf,3,2); % FA, #bands, f values
Mspgr = zeros(nf,3,2);

for IX = 1:2
    
    tmp = tissuepars;
    tmp.semi.f = f(IX);
    
    %%% Bieri-Scheffler correction
    tmp = Bieri_Scheffler_finite_correction(TBP,tau,TR,tmp);
    
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

%% Figure 3 from paper

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
    pl(1).LineStyle = ':';
    pl(1).Color = [0.2 0.4 1];
    hold on

    
    subplot(nr,nc,2+(jj-1)*2)
    plot(r2d(flips),abs(Mspgr(:,:,jj)),'linewidth',1.5)
    grid on
    xlabel('Flip angle, deg')
    ylabel('Signal / M_0')
    
    pl = get(gca,'Children');
    pl(3).Color = [0 0 0];
    pl(1).LineStyle = ':';
    pl(1).Color = [0.2 0.4 1];
    hold on
    
end

% add legend;
ll = legend('S_{1B}','S_{2B}','S_{3B}');
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
title('bSSFP','fontsize',16,'fontweight','bold')
text(-25,0.06,'f = 0','fontsize',18,'fontweight','bold','rotation',90,'fontangle','italic')

axes(gg(4))
title('SPGR','fontsize',16,'fontweight','bold')

axes(gg(3))
text(-25,0.06,'f = 1','fontsize',18,'fontweight','bold','rotation',90,'fontangle','italic')

%%% annotation
a1 = annotation('textarrow',[0.2786 0.2214],[0.6878 0.7715],'String','\DeltaMT','fontsize',13);
a2 = annotation('textarrow',[0.3014 0.2442],[0.2172 0.3009],'String','\DeltaihMT','fontsize',13);
a1.HeadStyle = 'cback1';
a2.HeadStyle = 'cback1';
a1.HeadWidth = 6;
a2.HeadWidth = 6;
%
setpospap([100 100 777 590])
print -dpng -r300 figs/simplefigure.png


%% Comparison for different frequency offsets and RMS B1, SSFP only

%%% Set up sequence

%%% generate some pulses
flipangle = d2r(30);

n=32;
b1rms = linspace(2,8,n);
deltaf= linspace(1000,16e3,n);

tau = 2e-3;
TR = 5e-3;
dt = 10e-6;

pulses = {};
b1sqrd = {};
b1pulse= {};
b1max=[];
for ii=1:n
    for jj=1:n
        % 3 bands
        [pulses{ii,jj,1},b1sqrd{ii,jj,1},b1pulse{ii,jj,1},TBP] = gen_MB_pulse(flipangle,tau,TR,b1rms(ii),deltaf(jj),'3','alpha',gauss_alpha,'dt',dt);
        % 2 bands
        [pulses{ii,jj,2},b1sqrd{ii,jj,2},b1pulse{ii,jj,2},TBP] = gen_MB_pulse(flipangle,tau,TR,b1rms(ii),deltaf(jj),'2+','alpha',gauss_alpha,'dt',dt);
    end
    b1max(ii,1)=max(abs(pulses{ii,n,1}));
    b1max(ii,2)=max(abs(pulses{ii,n,2})); %<-- should not depend on frequency offset
end

%%% Bieri-Scheffler correction
tmp = Bieri_Scheffler_finite_correction(TBP,tau,TR,tissuepars);

%%% Now simulate
ssfp_b1_delta = zeros([n n 3]);
for ii=1:n
    for jj=1:n
        
        dphi=0;
        
        %%% Symmetric (3-band) case
        ssfp_b1_delta(ii,jj,3) = ssSSFP_ihMT(flipangle,b1sqrd{ii,jj,1},deltaf(jj)*[-1 0 1],TR,tau,dphi,tmp);
        %%% Asymmetric (2-band) case
        ssfp_b1_delta(ii,jj,2) = ssSSFP_ihMT(flipangle,b1sqrd{ii,jj,2},deltaf(jj)*[-1 0 1],TR,tau,dphi,tmp);
        %%% Single band case
        ssfp_b1_delta(ii,jj,1) =ssSSFP_ihMT(flipangle,b1sqrd{ii,jj,1}.*[0 1 0],deltaf(jj)*[-1 0 1],TR,tau,dphi,tmp);
        
        disp([ii jj])
    end
end

%
figfp(1)
imagesc(deltaf*1e-3,b1rms,abs(ssfp_b1_delta(:,:,2)-ssfp_b1_delta(:,:,3))./abs(ssfp_b1_delta(:,:,3)),[0 0.15])
xlabel('\Delta/2\pi, kHz')
ylabel('B_1^{rms}, uT')
colorbar
set(gca,'fontsize',12)
title ('ihMT ratio vs pulse parameters')


%% Comparison for different frequency offsets and flip angles, SSFP only

%%% Set up sequence

%%% generate some pulses
flipangle = d2r(linspace(10,80,n));

n=32;
b1rmsfix = 5;
deltaf= linspace(1000,16e3,n);

tau = 2e-3;
TR = 5e-3;
dt = 10e-6;

pulses = {};
b1sqrd = {};
b1pulse= {};

for ii=1:n
    for jj=1:n
        % 3 bands
        [pulses{ii,jj,1},b1sqrd{ii,jj,1},b1pulse{ii,jj,1},TBP] = gen_MB_pulse(flipangle(ii),tau,TR,b1rmsfix,deltaf(jj),'3','alpha',gauss_alpha,'dt',dt);
        % 2 bands
        [pulses{ii,jj,2},b1sqrd{ii,jj,2},b1pulse{ii,jj,2},TBP] = gen_MB_pulse(flipangle(ii),tau,TR,b1rmsfix,deltaf(jj),'2+','alpha',gauss_alpha,'dt',dt);
    end
end

%%% Bieri-Scheffler correction
tmp = Bieri_Scheffler_finite_correction(TBP,tau,TR,tissuepars);

%%% Now simulate
ssfp_fa_delta = zeros([n n 3]);
for ii=1:n
    for jj=1:n
        
        dphi=0;
        
        %%% Symmetric (3-band) case
        ssfp_fa_delta(ii,jj,3) = ssSSFP_ihMT(flipangle(ii),b1sqrd{ii,jj,1},deltaf(jj)*[-1 0 1],TR,tau,dphi,tmp);
        %%% Asymmetric (2-band) case
        ssfp_fa_delta(ii,jj,2) = ssSSFP_ihMT(flipangle(ii),b1sqrd{ii,jj,2},deltaf(jj)*[-1 0 1],TR,tau,dphi,tmp);
        %%% Single band case
        ssfp_fa_delta(ii,jj,1) =ssSSFP_ihMT(flipangle(ii),b1sqrd{ii,jj,1}.*[0 1 0],deltaf(jj)*[-1 0 1],TR,tau,dphi,tmp);
        
        disp([ii jj])
    end
end

%%
figfp(2)
imagesc(deltaf*1e-3,r2d(flipangle),abs(ssfp_fa_delta(:,:,2)-ssfp_fa_delta(:,:,3))./abs(ssfp_fa_delta(:,:,1)))
xlabel('\Delta/2\pi, kHz')
ylabel('FA, deg')
colorbar
set(gca,'fontsize',12)
title ('ihMT ratio vs pulse parameters')


%% Design pulses for different TRs
tr = linspace(2e-3,10e-3,n);
b1max=zeros(n,n);
for jj=1:n
    for ii=1:n
        tmp=gen_MB_pulse(d2r(30),tau,tr(jj),b1rms(ii),8e3,'3','alpha',gauss_alpha,'dt',dt);
        b1max(ii,jj) = max(abs(tmp));
    end
end



%% Figure 7

nr=2;nc=3;
figfp(2)
subplot(nr,nc,2)
imagesc(deltaf*1e-3,b1rms,100*abs(ssfp_b1_delta(:,:,2)-ssfp_b1_delta(:,:,3)));
grid on
xlabel('\Delta/2\pi, kHz')
ylabel('B_1^{rms}, \muT')
colorbar
title ('\DeltaihMT / M_0')
text(17,8.5,'percent','fontweight','bold','fontsize',12,'rotation',0)
axis xy

subplot(nr,nc,1)
imagesc(deltaf*1e-3,flipangle*180/pi,100*abs(ssfp_fa_delta(:,:,2)-ssfp_fa_delta(:,:,3)));
grid on
xlabel('\Delta/2\pi, kHz')
ylabel('FA, deg')
colorbar
title ('\DeltaihMT')
text(17,85,'percent','fontweight','bold','fontsize',12,'rotation',0)
axis xy

subplot(nr,nc,5)
imagesc(deltaf*1e-3,b1rms,100*abs(ssfp_b1_delta(:,:,2)-ssfp_b1_delta(:,:,3))./abs(ssfp_b1_delta(:,:,1)));
grid on
xlabel('\Delta/2\pi, kHz')
ylabel('B_1^{rms}, \muT')
colorbar
%set(gca,'fontsize',12)
title ('ihMTR')
text(17,8.5,'percent','fontweight','bold','fontsize',12,'rotation',0)
axis xy

subplot(nr,nc,4)
% ihMT ratio
imagesc(deltaf*1e-3,flipangle*180/pi,100*abs(ssfp_fa_delta(:,:,2)-ssfp_fa_delta(:,:,3))./abs(ssfp_fa_delta(:,:,1)));
grid on

xlabel('\Delta/2\pi, kHz')
ylabel('FA, deg')
colorbar
%set(gca,'fontsize',12)
title ('ihMTR')
text(17,85,'percent','fontweight','bold','fontsize',12,'rotation',0)
axis xy

% B1 max
subplot(nr,nc,3)
[bb,trs]=meshgrid(tr*1e3,b1rms);
[cc,hh]=contourf(bb,trs,abs(b1max),5:5:40);
grid on
clabel(cc,hh,'Color',[1 1 1]);
hh.LineColor = [1 1 1];
hh.LineWidth = 2;
xlabel('TR, ms')
ylabel('B_1^{rms}, uT')
colorbar
% set(gca,'fontsize',12)
title ('Max B_1 (\muT)')
text(13,5.5,'\muT','fontweight','bold','fontsize',12,'rotation',0)

% add patch
pp1 = patch([2 2 4 4],[2 8 8 2],1);
pp1.FaceColor = [1 1 1];
pp1.EdgeColor = 'none';
pp1.FaceAlpha = 0.7;

pp2 = patch([4 4 10 10],[5 8 8 5],1);
pp2.FaceColor = [1 1 1];
pp2.EdgeColor = 'none';
pp2.FaceAlpha = 0.7;


% add text
text(-30,1,'(a)','fontweight','bold','fontsize',15)
text(-15,1,'(b)','fontweight','bold','fontsize',15)
text(-30,-7,'(c)','fontweight','bold','fontsize',15)
text(-15,-7,'(d)','fontweight','bold','fontsize',15)
text(0,1,'(e)','fontweight','bold','fontsize',15)


setpospap([100 254   1092 526])
print -dpng -r300 figs/whcontrastfig.png

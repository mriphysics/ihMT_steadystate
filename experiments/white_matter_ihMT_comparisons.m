%% 19-4-2019: This script performs simulations of bSSFP and SPGR for white matter
% Shaihan Malik, King's College London, 2019
% EDITS: 16-6-2019

%% Generates figures 4 and 8 from the paper

tissuepars = init_tissue('ic');


%%% Set up sequence

%%% generate some pulses
flips = deg2rad([1:2:30 35:5:80]);
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

%% Figure 4 from paper 

figfp(1)
nr=2;nc=2;

for jj=1:2
    
    subplot(nr,nc,1+(jj-1)*nc)
    plot(rad2deg(flips),abs(Mssfp(:,:,jj)),'linewidth',1.5);
    grid on
    xlabel('Flip angle, deg')
    ylabel('Signal / M_0')
    pl = get(gca,'Children');
    pl(3).Color = [0 0 0];
    pl(1).LineStyle = '--';
    pl(1).Color = [0 0 1];
    pl(2).Color = [0 1 0];
    
    hold on

    
    subplot(nr,nc,2+(jj-1)*nc)
    plot(rad2deg(flips),abs(Mspgr(:,:,jj)),'linewidth',1.5)
    grid on
    xlabel('Flip angle, deg')
    ylabel('Signal / M_0')
    
    pl = get(gca,'Children');
    pl(3).Color = [0 0 0];
    pl(1).LineStyle = '--';
    pl(1).Color = [0 0 1];
    pl(2).Color = [0 1 0];
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
x2 = [rad2deg(flips) flip(rad2deg(flips),2)];
fh=fill(x2, inBetween, 'b');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.2;

axes(gg(4))
inBetween = [abs(Mspgr(:,1,1)); flip(abs(Mspgr(:,2,1)),1)];
x2 = [rad2deg(flips) flip(rad2deg(flips),2)];
fh=fill(x2, inBetween, 'b');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.2;

axes(gg(3))
inBetween = [abs(Mssfp(:,2,2)); flip(abs(Mssfp(:,3,2)),1)];
x2 = [rad2deg(flips) flip(rad2deg(flips),2)];
fh=fill(x2, inBetween, 'g');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.2;

axes(gg(2))
inBetween = [abs(Mspgr(:,2,2)); flip(abs(Mspgr(:,3,2)),1)];
x2 = [rad2deg(flips) flip(rad2deg(flips),2)];
fh=fill(x2, inBetween, 'g');
fh.EdgeColor = 'none';
fh.FaceAlpha=0.2;

%%% text labels
axes(gg(5))
title('bSSFP','fontsize',16,'fontweight','bold')
text(-25,0.06,'\delta = 0','fontsize',18,'fontweight','bold','rotation',90,'fontangle','italic')

axes(gg(4))
title('SPGR','fontsize',16,'fontweight','bold')

axes(gg(3))
text(-25,0.06,'\delta = 1','fontsize',18,'fontweight','bold','rotation',90,'fontangle','italic')

%%% annotation
a1 = annotation('textarrow',[0.2786 0.2214],[0.6878 0.7715],'String','\DeltaMT','fontsize',13);
a2 = annotation('textarrow',[0.3181 0.2609],[0.2121 0.2958],'String','\DeltaihMT','fontsize',13);
a1.HeadStyle = 'cback1';
a2.HeadStyle = 'cback1';
a1.HeadWidth = 6;
a2.HeadWidth = 6;
%
for ii=2:5
    gg(ii).Color = [1 1 1]*0.95;
end
setpospap([100 100 777 590])
set(gcf, 'InvertHardcopy', 'off')
print -dpng -r300 figs/simplefigure.png

%% Add supporting information figure on ihMTR and B1rms
figfp(2)

subplot(3,1,1)
plot(rad2deg(flips),100*(abs(Mspgr(:,2,2))-abs(Mspgr(:,3,2))),'linewidth',1.5)
grid on
hold on
plot(rad2deg(flips),100*(abs(Mssfp(:,2,2))-abs(Mssfp(:,3,2))),'linewidth',1.5)
ylim([0 1.4])
xlabel('Flip angle, deg')
ylabel('\DeltaihMT (percent of M_0)')
title('\DeltaihMT','fontsize',14,'fontweight','bold')
ll=legend('SPGR','bSSFP');ll.FontSize=11;

subplot(3,1,2)
plot(rad2deg(flips),100*(abs(Mspgr(:,2,2))-abs(Mspgr(:,3,2)))./abs(Mspgr(:,1,2)),'linewidth',1.5)
grid on
hold on
plot(rad2deg(flips),100*(abs(Mssfp(:,2,2))-abs(Mssfp(:,3,2)))./abs(Mssfp(:,1,2)),'linewidth',1.5)

xlabel('Flip angle, deg')
ylabel('ihMTR  (% of S_{1B})')
title('ihMTR','fontsize',14,'fontweight','bold')
ll=legend('SPGR','bSSFP');ll.FontSize=11;



subplot(3,1,3)
b1offres=[];for ii=1:nf,b1offres(ii)=sqrt(b1sqrd{ii,2}(3)*tau/TR);end
plot(rad2deg(flips),b1offres,'linewidth',1.5)
grid on
xlabel('Flip angle, deg')
ylabel('$\sqrt<B_1^2(\Delta\neq 0)>,\mu T$','interpreter','Latex')
title('Off-resonance B_1^{rms}','fontsize',14,'fontweight','bold')
ylim([0 6]);

setpospap([00 00 400 700])
print -dpng -r300 figs/simplefigure_supportinginfo.png

%% Comparison for different frequency offsets and RMS B1, SSFP only

%%% Set up sequence

%%% generate some pulses
flipangle = deg2rad(30);

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
flipangle = deg2rad(linspace(10,80,n));

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
imagesc(deltaf*1e-3,rad2deg(flipangle),abs(ssfp_fa_delta(:,:,2)-ssfp_fa_delta(:,:,3))./abs(ssfp_fa_delta(:,:,1)))
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
        tmp=gen_MB_pulse(deg2rad(30),tau,tr(jj),b1rms(ii),8e3,'3','alpha',gauss_alpha,'dt',dt);
        b1max(ii,jj) = max(abs(tmp));
    end
end



%% Figure 8

nr=2;nc=3;
fs=13;
figfp(2)
subplot(nr,nc,2)
imagesc(deltaf*1e-3,b1rms,100*abs(ssfp_b1_delta(:,:,2)-ssfp_b1_delta(:,:,3)));
grid on
xlabel('\Delta/2\pi, kHz')
ylabel('B_1^{rms}, \muT')
colorbar
title ('\DeltaihMT','fontsize',fs)
text(21.5,4,'% of M_0','fontweight','bold','fontsize',12,'rotation',90)
axis xy

subplot(nr,nc,1)
imagesc(deltaf*1e-3,flipangle*180/pi,100*abs(ssfp_fa_delta(:,:,2)-ssfp_fa_delta(:,:,3)));
grid on
xlabel('\Delta/2\pi, kHz')
ylabel('FA, deg')
colorbar
title ('\DeltaihMT','fontsize',fs)
text(21.5,35,'% of M_0','fontweight','bold','fontsize',12,'rotation',90)
axis xy

subplot(nr,nc,5)
imagesc(deltaf*1e-3,b1rms,100*abs(ssfp_b1_delta(:,:,2)-ssfp_b1_delta(:,:,3))./abs(ssfp_b1_delta(:,:,1)));
grid on
xlabel('\Delta/2\pi, kHz')
ylabel('B_1^{rms}, \muT')
colorbar
%set(gca,'fontsize',12)
title ('ihMTR','fontsize',fs)
text(21.5,4,'% of S_{1B}','fontweight','bold','fontsize',12,'rotation',90)
axis xy

subplot(nr,nc,4)
% ihMT ratio
imagesc(deltaf*1e-3,flipangle*180/pi,100*abs(ssfp_fa_delta(:,:,2)-ssfp_fa_delta(:,:,3))./abs(ssfp_fa_delta(:,:,1)));
grid on

xlabel('\Delta/2\pi, kHz')
ylabel('FA, deg')
colorbar
%set(gca,'fontsize',12)
title ('ihMTR','fontsize',fs)
text(21.5,35,'% of S_{1B}','fontweight','bold','fontsize',12,'rotation',90)
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
title ('Max B_1 (\muT)','fontsize',fs)
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

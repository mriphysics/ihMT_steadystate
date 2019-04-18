%%% Test the SSFP simulations

tissuepars = init_tissue();


%%% Set up sequence

%%% generate some pulses
flips = d2r(5:5:70);
nf = length(flips);
tau = 2e-3;
TR = 5e-3;
b1_rms = 4;
delta = 7e3;
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


%% bSSFP steady state

%%% Extra pars needed for SSFP
dphi=0;


% loop over the different flips
M = zeros(nf,3);
Mz = zeros(nf,2,4);
Delta_Hz = [-delta 0 delta];
for ii=1:nf
    %%% Symmetric (3-band) case
    [M(ii,1),Mz(ii,1,:)] = ssSSFP_ihMT(flips(ii),b1sqrd{ii,1},Delta_Hz,TR,tau,dphi,tissuepars);
    %%% Asymmetric (2-band) case
    [M(ii,2),Mz(ii,2,:)] = ssSSFP_ihMT(flips(ii),b1sqrd{ii,2},Delta_Hz,TR,tau,dphi,tissuepars);
    %%% Single band case
    M(ii,3) =ssSSFP_ihMT(flips(ii),b1sqrd{ii,1}.*[0 1 0],Delta_Hz,TR,tau,dphi,tissuepars);
end

figfp(21);
subplot(3,1,1)
plot(r2d(flips),abs(M),'linewidth',2)
legend('3-band','2-band','1-band')
grid on
set(gca,'fontsize',12)
xlabel('flip')
ylabel('signal')
title('Predicted signals')
% ylim([0 0.08]);

subplot(3,1,2)
plot(r2d(flips),abs(M(:,2)-M(:,1)),'linewidth',2)
legend('CNR')
grid on
set(gca,'fontsize',12)
xlabel('flip')
ylabel('signal')
title('ihMT difference')

subplot(3,1,3)
cmap=colormap(lines);
hold on
ll={'--','-'};
for ii=1:2
    for jj=1:3
        plot(r2d(flips),squeeze(Mz(:,ii,jj)),ll{ii},'linewidth',1.2,'color',cmap(jj,:));
    end
end
legend('M_{free}^{Zeeman} (3-band)','M_{bound}^{Zeeman} (3-band)','M_{bound}^{Dipolar} (3-band)');%,...
    %'M_{free}^{Zeeman} (2-band)','M_{bound}^{Zeeman} (2-band)','M_{bound}^{Dipolar} (2-band)')

grid on
set(gca,'fontsize',12)
xlabel('flip')
ylabel('signal')
title('Predicted Mz terms')


%% Full integration

%%% Extra pars needed for SSFP
dphi=0;


% loop over the different flips
Mint = zeros(nf,3);
Delta_Hz = [-delta 0 delta];
for ii=1:nf
    %%% Symmetric (3-band) case
    [Mint(ii,1)] = ssSSFP_ihMT_integrate(b1pulse{ii,1},dt,Delta_Hz,TR,dphi,tissuepars);
    %%% Asymmetric (2-band) case
    [Mint(ii,2)] = ssSSFP_ihMT_integrate(b1pulse{ii,2},dt,Delta_Hz,TR,dphi,tissuepars);
    %%% Single band case
    Mint(ii,3) =ssSSFP_ihMT_integrate(b1pulse{ii,1}.*[0 1 0],dt,Delta_Hz,TR,dphi,tissuepars);
end

figfp(22);
plot(r2d(flips),abs(M),'linewidth',2)
hold on
plot(r2d(flips),abs(Mint),'--','linewidth',2)

legend('3-band','2-band','1-band','INT 3-band','INT 2-band','INT 1-band')
grid on
set(gca,'fontsize',12)
xlabel('flip')
ylabel('signal')
title('Predicted signals')
% ylim([0 0.08]);

%% Steady-state with Bieri-Scheffler correction
%%% Extra pars needed for SSFP
dphi=0;

% compute new R2
TBP=2; % from power spectrum
tissuepars_mod = Bieri_Scheffler_finite_correction(TBP,tau,TR,tissuepars);


% loop over the different flips
Mhat = zeros(nf,3);
Delta_Hz = [-delta 0 delta];
for ii=1:nf
    %%% Symmetric (3-band) case
    [Mhat(ii,1)] = ssSSFP_ihMT(flips(ii),b1sqrd{ii,1},Delta_Hz,TR,tau,dphi,tissuepars_mod);
    %%% Asymmetric (2-band) case
    [Mhat(ii,2)] = ssSSFP_ihMT(flips(ii),b1sqrd{ii,2},Delta_Hz,TR,tau,dphi,tissuepars_mod);
    %%% Single band case
    Mhat(ii,3) =ssSSFP_ihMT(flips(ii),b1sqrd{ii,1}.*[0 1 0],Delta_Hz,TR,tau,dphi,tissuepars_mod);
end


figfp(23);
for ii=1:3
    subplot(1,3,ii)
    plot(r2d(flips),abs(M(:,ii)))
    hold on
    plot(r2d(flips),abs(Mint(:,ii)))
    plot(r2d(flips),abs(Mhat(:,ii)))
    
    legend('SS','int','SS corrected')
    grid on
    set(gca,'fontsize',12)
    xlabel('flip')
    ylabel('signal')
    title('Predicted signals')
end


%% Now consider a fixed flip angle, look at different pulse duration and dipolar relaxation time

%%% generate some pulses
flipangle = d2r(60);

% range of tau and dipolar T1
n=32;
tau = linspace(0.2e-3,2.5e-3,n);
T1D = linspace(0.1e-3,30e-3,n);
T1D = logspace(-4,-2,n);

TR = 5e-3;
b1_rms = 5;
delta = 7e3;
dt = 10e-6;
TBP = 2;

pulses = {};
b1sqrd = {};
b1pulse= {};
for ii=1:n
    % 3 bands
    [pulses{ii,1},b1sqrd{ii,1},b1pulse{ii,1}] = gen_MB_pulse(flipangle,tau(ii),TR,b1_rms,delta,'3','sigma',2,'dt',dt); 
    % 2 bands
    [pulses{ii,2},b1sqrd{ii,2},b1pulse{ii,2}] = gen_MB_pulse(flipangle,tau(ii),TR,b1_rms,delta,'2+','sigma',2,'dt',dt);
end


%%% Now compute a grid of points
M = zeros([n n 3 3]); % T1D, tau, method, number of bands

for ii=1:n

    % generate new tissue parameter set
    tmp = init_tissue();
    tmp.semi.R1D = 1/T1D(ii);

    % Now loop over pulses (i.e. different tau)
    for jj=1:n
        
        %==== Standard Steady state method ====
        %%% Symmetric (3-band) case
        M(ii,jj,1,1) = ssSSFP_ihMT(flipangle,b1sqrd{jj,1},Delta_Hz,TR,tau(jj),dphi,tmp);
        %%% Asymmetric (2-band) case
        M(ii,jj,1,2) = ssSSFP_ihMT(flipangle,b1sqrd{jj,2},Delta_Hz,TR,tau(jj),dphi,tmp);
        %%% Single band case
        M(ii,jj,1,3) =ssSSFP_ihMT(flipangle,b1sqrd{jj,1}.*[0 1 0],Delta_Hz,TR,tau(jj),dphi,tmp);
        
        %==== Full integration ====
        %%% Symmetric (3-band) case
        M(ii,jj,2,1) = ssSSFP_ihMT_integrate(b1pulse{jj,1},dt,Delta_Hz,TR,dphi,tmp);
        %%% Asymmetric (2-band) case
        M(ii,jj,2,2) = ssSSFP_ihMT_integrate(b1pulse{jj,2},dt,Delta_Hz,TR,dphi,tmp);
        %%% Single band case
        M(ii,jj,2,3) =ssSSFP_ihMT_integrate(b1pulse{jj,1}.*[0 1 0],dt,Delta_Hz,TR,dphi,tmp);
        
        %=== Bieri Scheffler Correction
        tmpC = Bieri_Scheffler_finite_correction(TBP,tau(jj),TR,tmp);

        %%% Symmetric (3-band) case
        M(ii,jj,3,1) = ssSSFP_ihMT(flipangle,b1sqrd{jj,1},Delta_Hz,TR,tau(jj),dphi,tmpC);
        %%% Asymmetric (2-band) case
        M(ii,jj,3,2) = ssSSFP_ihMT(flipangle,b1sqrd{jj,2},Delta_Hz,TR,tau(jj),dphi,tmpC);
        %%% Single band case
        M(ii,jj,3,3) =ssSSFP_ihMT(flipangle,b1sqrd{jj,1}.*[0 1 0],Delta_Hz,TR,tau(jj),dphi,tmpC);
        
       disp([ii jj]) 
    end
end


%%
figfp(1)
nb = [3 2 1];
for ii=1:3
    subplot(2,3,ii)
%     imsjm((abs(M(:,:,2,ii))-abs(M(:,:,1,ii)))./abs(M(:,:,2,ii)),[0 0.1])
    [cc,hh]=contourf(1e3*tau,T1D,100*abs(abs(M(:,:,2,nb(ii)))-abs(M(:,:,1,nb(ii))))./abs(M(:,:,2,nb(ii))),0:1:10);
    title(sprintf('%d band, error',ii))
    set(gca,'Yscale','log')
    grid on
    caxis([0 10])
    hh.LineColor = [1 1 1];
    
    xlabel('Pulse duration \tau (ms)')
    ylabel('Dipolar relaxation time T_{1D}^s (s)')
    
    subplot(2,3,ii+3)
    %imsjm(abs(abs(M(:,:,2,ii))-abs(M(:,:,3,ii)))./abs(M(:,:,2,ii)),[-0.005 0.1])
    [cc,hh]=contourf(1e3*tau,T1D,100*abs(abs(M(:,:,2,nb(ii)))-abs(M(:,:,3,nb(ii))))./abs(M(:,:,2,nb(ii))),0:1:10);
    title(sprintf('%d band, error with correction',ii))
    set(gca,'Yscale','log')
    caxis([0 10])
    xlabel('Pulse duration \tau (ms)')
    ylabel('Dipolar relaxation time T_{1D}^s (s)')
    hh.LineColor = [1 1 1];
    grid on
       

    colormap hot
end

cc = colorbar;
cc.Position = [0.93 0.3 0.02 0.5];


%% Off-resonance properties


tissuepars = init_tissue();

flipngle = d2r(50);
tau = 2e-3;
TR = 5e-3;
b1_rms = 4;
delta = 7e3;
dt = 10e-6;


% 3 bands
[~,b1sqrd3] = gen_MB_pulse(flipangle,tau,TR,b1_rms,delta,'3','sigma',2,'dt',dt);
% 2 bands
[~,b1sqrd2] = gen_MB_pulse(flipangle,tau,TR,b1_rms,delta,'2+','sigma',2,'dt',dt);


n=64;
dphi = linspace(-pi*1.5,pi*1.5,n);
Mphi= zeros([n 3]);

for ii=1:n
    %%% Symmetric (3-band) case
    Mphi(ii,1) = ssSSFP_ihMT(flips(ii),b1sqrd3,Delta_Hz,TR,tau,dphi(ii),tissuepars);
    %%% Asymmetric (2-band) case
    Mphi(ii,2) = ssSSFP_ihMT(flips(ii),b1sqrd2,Delta_Hz,TR,tau,dphi(ii),tissuepars);
    %%% Single band case
    Mphi(ii,3) =ssSSFP_ihMT(flips(ii),b1sqrd3.*[0 1 0],Delta_Hz,TR,tau,dphi(ii),tissuepars);
end

figfp(1)
plot(dphi,abs(Mphi))

%% Loop over different B1rms

flipangle = d2r(30);

tau = 2e-3;
TR = 5e-3;
b1_rms = 1:0.2:20;
nf = length(b1_rms);

delta = 7e3;
dt = 10e-6;

pulses = {};
b1sqrd = {};
b1pulse= {};
for ii=1:nf
    % 3 bands
    [pulses{ii,1},b1sqrd{ii,1},b1pulse{ii,1}] = gen_MB_pulse(flipangle,tau,TR,b1_rms(ii),delta,'3','sigma',2,'dt',dt); 
    % 2 bands
    [pulses{ii,2},b1sqrd{ii,2},b1pulse{ii,2}] = gen_MB_pulse(flipangle,tau,TR,b1_rms(ii),delta,'2+','sigma',2,'dt',dt);
end



%%% Extra pars needed for SSFP
dphi=0;


% loop over the different flips
M = zeros(nf,3);
Mz = zeros(nf,2,4);
Delta_Hz = [-delta 0 delta];
for ii=1:nf
    %%% Symmetric (3-band) case
    [M(ii,1),Mz(ii,1,:)] = ssSSFP_ihMT(flipangle,b1sqrd{ii,1},Delta_Hz,TR,tau,dphi,tissuepars);
    %%% Asymmetric (2-band) case
    [M(ii,2),Mz(ii,2,:)] = ssSSFP_ihMT(flipangle,b1sqrd{ii,2},Delta_Hz,TR,tau,dphi,tissuepars);
    %%% Single band case
    M(ii,3) =ssSSFP_ihMT(flipangle,b1sqrd{ii,1}.*[0 1 0],Delta_Hz,TR,tau,dphi,tissuepars);
end

figfp(21);
subplot(2,1,1)
plot(b1_rms,abs(M),'linewidth',2)
legend('3-band','2-band','1-band')
grid on
set(gca,'fontsize',12)
xlabel('B_1^{rms}')
ylabel('signal')
title('Predicted signals')
% ylim([0 0.08]);

subplot(2,1,2)
plot(b1_rms,abs(M(:,2)-M(:,1)),'linewidth',2)
legend('CNR')
grid on
set(gca,'fontsize',12)
xlabel('B_1^{rms}')
ylabel('signal')
title('ihMT difference')

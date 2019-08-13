%%% Script to re-generate all fitting of the experimental phantom data. The
%%% phantom consisted of water, BSA, PL-161, and hair conditioner.


%% Define SSFP pulse sequence
seq_pars = struct;
seq_pars.flips = deg2rad(10:10:80);
nf = length(seq_pars.flips);
seq_pars.tau = 2.2e-3;
seq_pars.TR = 5e-3;
seq_pars.b1_rms = 4.15;
seq_pars.delta = 8e3;
dt = 10e-6;

seq_pars.b1sqrd = {};

for ii=1:nf
    % 3 bands
    [~,seq_pars.b1sqrd{ii,3},~,seq_pars.TBP] = gen_MB_pulse(seq_pars.flips(ii),seq_pars.tau,seq_pars.TR,seq_pars.b1_rms,seq_pars.delta,'3','alpha',3,'dt',dt); 
    % 2 bands
    [~,seq_pars.b1sqrd{ii,2},~] = gen_MB_pulse(seq_pars.flips(ii),seq_pars.tau,seq_pars.TR,seq_pars.b1_rms,seq_pars.delta,'2+','alpha',3,'dt',dt);
    % 1 bands
    seq_pars.b1sqrd{ii,1} = seq_pars.b1sqrd{ii,2}.*[0 1 0]; 
end
seq_pars.Delta_Hz = seq_pars.delta * [-1 0 1];
seq_pars.dphi=0;

%% Define SPGR pulse sequence
seq_pars_spgr = seq_pars;
seq_pars_spgr.flips = deg2rad(2:2:16);
nf = length(seq_pars_spgr.flips);

seq_pars_spgr.b1sqrd = {};

for ii=1:nf
    % 3 bands
    [~,seq_pars_spgr.b1sqrd{ii,3},~,TBP] = gen_MB_pulse(seq_pars_spgr.flips(ii),seq_pars_spgr.tau,seq_pars_spgr.TR,seq_pars_spgr.b1_rms,seq_pars_spgr.delta,'3','alpha',3,'dt',dt); 
    % 2 bands
    [~,seq_pars_spgr.b1sqrd{ii,2},~] = gen_MB_pulse(seq_pars_spgr.flips(ii),seq_pars_spgr.tau,seq_pars_spgr.TR,seq_pars_spgr.b1_rms,seq_pars_spgr.delta,'2+','alpha',3,'dt',dt);
    % 1 bands
    seq_pars_spgr.b1sqrd{ii,1} = seq_pars_spgr.b1sqrd{ii,2}.*[0 1 0]; 
end
seq_pars_spgr.Delta_Hz = seq_pars.delta * [-1 0 1];

% add echo time
seq_pars_spgr.TE = 2.5e-3;

%% Load in data to fit

load bin/fitdata.mat
IX = 3; %PL161
IX = 2; % BSA
IX = 4; % Hair cond
IX = 5; % MnCl2

IX = [5 2 3 4];

%%% Get the data for each phantom into a cell array
xdata = {};
sdata = {};
for jj=1:4
    tmp = squeeze(kdata{IX(jj),1});
    % average the 2+ and 2- data
    tmp(:,:,2) = 0.5*(tmp(:,:,2)+tmp(:,:,3));
    tmp(:,:,3) = [];
    
    %%% Now the SPGR data
    tmp2 = squeeze(kdata{IX(jj),2});
    
    % average the 2+ and 2- data
    tmp2(:,:,2) = 0.5*(tmp2(:,:,2)+tmp2(:,:,3));
    tmp2(:,:,3) = [];
    
    % reorder and scale down by 1000
    xdata{jj,1} = permute(tmp/1000,[2 3 1]);
    xdata{jj,2} = permute(tmp2/1000,[2 3 1]);
 
    % standard deviation
    sdata{jj,1} = std(xdata{jj,1},1,3);
    sdata{jj,2} = std(xdata{jj,2},1,3);
    % also store number of samples
    sdata{jj,3} = size(xdata{jj,1},3);
    
    % average over the whole phantom 
    xdata{jj,1} = mean(xdata{jj,1},3);
    xdata{jj,2} = mean(xdata{jj,2},3);
    
end


%% Combined Cost Function and fitting. 


%%% Fix bounds for each one
%       [R1f R2 M0f  R1sz  R1D  f  T2s(us)  k   sf    b1sf]
lb={};
lb{1} = [0.2 1  0    0.8   30   0  10       10  0   1.0]; %<--- MnCl2
lb{2} = [0.2 1  0    0.8   30   0  10       10  0   1.0]; %<--- BSA
lb{3} = [0.2 1  0    0.8   30   1  10       10  0   1.0]; %<--- PL161
lb{4} = [0.2 1  0    0.8   30   0  10       10  0   1.0]; %<--- HC
ub={};
ub{1} = [3   40 0    0.8   30   0  10       10  10   1.0];%<-- just fit free pool, water
ub{2} = [3   40 1    50    66   0  25       100 10   1.0]; %<--- BSA
ub{3} = [3   40 1    50    66   1  25       100 10   1.0]; %<--- PL161
ub{4} = [3   40 1    50    66   1  25       100 10   1.0]; %<--- HC



%%% These are previously found best fits, use as starting points here just
%%% for demonstration
xs_HC = [0.5955    5.9999    0.0675    2.9780   42.8042    0.7139   25.8479   67.8720    5.1397    1.0000];
xs_BSA= [0.5635   18.8122    0.0774    7.7545   47.8069         0   15.2003   70.6874    4.8311    1.0000];
xs_PL = [0.5515   13.6859    0.1516    4.4111   48.2443    1.0000   17.8722   46.7160    4.9972    1.0000];%<-- fix f=1
xs_Mn = [0.5205    1.6267         0    4.8914   48.2782    0.4920   17.8562   46.6711    5.2583    1.0000];

%%% Use above best solutions as starting points for fit
x0 = {};
x0{1} = xs_Mn;
x0{2} = xs_BSA;
x0{3} = xs_PL;
x0{4} = xs_HC;

%%% Single forward model function
flag = 1; %<-- Gaussian line
fwd_model = @(x)(cat(2,ssfp_ihmt_fit_fwd_model(x,seq_pars,flag),spgr_ihmt_fit_fwd_model(x,seq_pars_spgr,flag)));


%%% Now loop over phantoms

for ix = 1:4
    
    % data
    data0 = cat(2,xdata{ix,1},xdata{ix,2});
    
    % cost function given this data
    costf = @(x)(norm(fwd_model(x)-data0)^2);
    
    % cost function in percentage terms
    costf_prcnt = @(x)(100*sqrt(costf(x))/norm(data0));
    
    options = optimoptions('fmincon');
    options = optimoptions(options,'Display', 'iter');
    options = optimoptions(options,'PlotFcns', {  @optimplotx @optimplotfunccount @optimplotfval });
    
    [xs,xf,xflg,output] = fmincon(costf,x0{ix},[],[],[],[],lb{ix},ub{ix},[],options);
    
    % store solution
    sol{ix} = xs;
    
    fprintf('Phantom %d nrmse %1.2f \n',ix,costf_prcnt(xs))
    nrmse(ix) = costf_prcnt(xs);
    
end


%% Generate figure with all of the fits

titles = {'Water with MnCl_2','BSA','PL-161','Hair Conditioner'};
nr=2;nc=2;
ms=4;

figfp(1)

for jj=1:4
    
    %%% Get the data again, this time don't average S2B+ and S2B- or
    %%% rescale by 1000
    xd1 = squeeze(mean(kdata{IX(jj),1},1));
    xd2 = squeeze(mean(kdata{IX(jj),2},1));
      
    subplot(nr,nc,jj)
    
    %%% Plot Data
    tmp1=plot(seq_pars.flips*180/pi,xd1,'o','markersize',ms);
    for ii=1:4
        tmp1(ii).MarkerFaceColor = tmp1(ii).Color;
    end
    hold on
    tmp2=plot(seq_pars_spgr.flips*180/pi,xd2,'V','markersize',ms);
    for ii=1:4
        tmp2(ii).MarkerFaceColor = tmp1(ii).Color;
        tmp2(ii).Color = tmp1(ii).Color;
    end
    grid on
    xlabel('Flip angle, degrees')
    ylabel('Signal (au)')
    
    if jj==2
        ll=legend('bSSFP: 1B','bSSFP: 2^+B','bSSFP: 2^-B','bSSFP: 3B','SPGR:  1B','SPGR:  2^+B','SPGR:  2^-B','SPGR:  3B','AutoUpdate','off');
        ll.Location = 'EastOutside';
        ll.FontSize = 9;
    end
    
    %%% Plot fits, scale up by 1000 to match unscaled data
    plot(seq_pars.flips*180/pi,1000*ssfp_ihmt_fit_fwd_model(sol{jj},seq_pars,flag),'-','Color',[0 0 0])
    plot(seq_pars_spgr.flips*180/pi,1000*spgr_ihmt_fit_fwd_model(sol{jj},seq_pars_spgr,flag),'-','Color',[0 0 0])
    
    %%% Label NRMSE
    yl=get(gca,'YLim');
    ylim([0 1.1*yl(2)]);
    text(40,yl(1)+(yl(2)-yl(1))/3.8,sprintf('NRMSE = %1.2f %%',nrmse(jj)),'fontweight','bold','fontsize',12)
    
    title(titles{jj},'fontsize',12)
end
setpospap([100 100 850 550])

gg = get(gcf,'Children');

gg(5).Position = [0.07 0.58 0.32 0.34];
gg(2).Position = [0.07 0.11 0.32 0.34];
gg(4).Position = [0.49 0.58 0.32 0.34];
gg(1).Position = [0.49 0.11 0.32 0.34];
gg(3).Position = [0.8296    0.4095    0.1559    0.2314];


% print -dpng -r300 figs/ihMT_phantom_plots.png

%% Alternative version with error bars

titles = {'Water with MnCl_2','BSA','PL-161','Hair Conditioner'};
nr=2;nc=2;
ms=3;

figfp(1)

for jj=1:4
    
    %%% Get the data again, this time don't average S2B+ and S2B- or
    %%% rescale by 1000
    xd1 = squeeze(mean(kdata{IX(jj),1},1));
    xd2 = squeeze(mean(kdata{IX(jj),2},1));

    %%% standard error
    n = size(kdata{IX(jj),1},1);
    CI95 = tinv([0.025 0.975],n-1);
    sd1 = CI95(2)*squeeze(std(kdata{IX(jj),1},[],1))/sqrt(n);
    sd2 = CI95(2)*squeeze(std(kdata{IX(jj),2},[],1))/sqrt(n);
    
    subplot(nr,nc,jj)
    
    %%% Plot Data
    %tmp1=plot(seq_pars.flips*180/pi,xd1,'o','markersize',ms);
    tmp1 = errorbar(repmat(seq_pars.flips'*180/pi,1,4),xd1,sd1,'s','markersize',ms);
    for ii=1:4
        tmp1(ii).MarkerFaceColor = tmp1(ii).Color;
    end
    hold on
    %tmp2=plot(seq_pars_spgr.flips*180/pi,xd2,'V','markersize',ms);
    tmp2 = errorbar(repmat(seq_pars_spgr.flips'*180/pi,1,4),xd2,sd2,'^','markersize',ms);
    for ii=1:4
        tmp2(ii).MarkerFaceColor = tmp1(ii).Color;
        tmp2(ii).Color = tmp1(ii).Color;
    end
    grid on
    xlabel('Flip angle, degrees')
    ylabel('Signal (au)')
    
    if jj==2
        ll=legend('bSSFP: 1B','bSSFP: 2^+B','bSSFP: 2^-B','bSSFP: 3B','SPGR:  1B','SPGR:  2^+B','SPGR:  2^-B','SPGR:  3B','AutoUpdate','off');
        ll.Location = 'EastOutside';
        ll.FontSize = 9;
    end
    
    %%% Plot fits, scale up by 1000 to match unscaled data
    plot(seq_pars.flips*180/pi,1000*ssfp_ihmt_fit_fwd_model(sol{jj},seq_pars,flag),'-','Color',[0 0 0])
    plot(seq_pars_spgr.flips*180/pi,1000*spgr_ihmt_fit_fwd_model(sol{jj},seq_pars_spgr,flag),'-','Color',[0 0 0])
    
    %%% Label NRMSE
    yl=get(gca,'YLim');
    ylim([0 1.1*yl(2)]);
    text(40,yl(1)+(yl(2)-yl(1))/3.8,sprintf('NRMSE = %1.2f %%',nrmse(jj)),'fontweight','bold','fontsize',12)
    
    title(titles{jj},'fontsize',12)
end
setpospap([100 100 850 550])

gg = get(gcf,'Children');

gg(5).Position = [0.07 0.58 0.32 0.34];
gg(2).Position = [0.07 0.11 0.32 0.34];
gg(4).Position = [0.49 0.58 0.32 0.34];
gg(1).Position = [0.49 0.11 0.32 0.34];
gg(3).Position = [0.8296    0.4095    0.1559    0.2314];





%% Add a plot of the ihMT ratio and ihMT difference for measured and simulated data based on 
% fitted parameters

titles = {'Water with MnCl_2','BSA','PL-161','Hair Conditioner'};
nr=2;nc=2;
ms=4;

figfp(1)

for jj=1:4
   
    subplot(nr,nc,jj)
    
    %%% Plot Data
    hold on
    tmp3=plot(seq_pars.flips*180/pi,100*(xdata{jj,1}(:,1)-xdata{jj,1}(:,2))./xdata{jj,1}(:,1),'o','markersize',ms,'color',[0 0 1]);
    tmp3.MarkerFaceColor = [0 0 1];
    tmp4=plot(seq_pars_spgr.flips*180/pi,100*(xdata{jj,2}(:,1)-xdata{jj,2}(:,2))./xdata{jj,2}(:,1),'^','markersize',ms,'color',[0 0 1]);
    tmp4.MarkerFaceColor = [0 0 1];
    tmp1=plot(seq_pars.flips*180/pi,100*(xdata{jj,1}(:,2)-xdata{jj,1}(:,3))./xdata{jj,1}(:,1),'o','markersize',ms,'color',[0 1 0]);
    tmp1.MarkerFaceColor = [0 1 0];
    tmp2=plot(seq_pars_spgr.flips*180/pi,100*(xdata{jj,2}(:,2)-xdata{jj,2}(:,3))./xdata{jj,2}(:,1),'^','markersize',ms,'color',[0 1 0]);
    tmp2.MarkerFaceColor = [0 1 0];
    
    grid on
    grid on
    xlabel('Flip angle, degrees')
    ylabel('Ratio (% of S_{1B})')
    
    if jj==1
        ll=legend('MTR, bSSFP','MTR, SPGR','ihMTR, bSSFP','ihMTR, SPGR','AutoUpdate','off');
        %ll.Location = 'EastOutside';
        ll.FontSize = 12;
    end
    yl = get(gca,'YLim');
    ylim([yl(1) 60])
    
    gg=gca;
    gg.XAxisLocation = 'origin';

    title(titles{jj},'fontsize',12)
end
setpospap([100 100 850 550])


%% Repeat above but with proper error propagation for MTR and ihMTR

titles = {'Water with MnCl_2','BSA','PL-161','Hair Conditioner'};
nr=2;nc=2;
ms=4;

figfp(1)

for jj=1:4
   
    subplot(nr,nc,jj)
    
    hold on 
    grid on
    
    for kk=1:2 %<--- loop over sequence type 1=bSSFP, 2=SPGR
        
        ihMTR = 100*(xdata{jj,kk}(:,2)-xdata{jj,kk}(:,3))./xdata{jj,kk}(:,1);
        var_ihMTR = ihMTR.^2 .* ( (sdata{jj,kk}(:,2).^2+sdata{jj,kk}(:,3).^2)./(xdata{jj,kk}(:,2)-xdata{jj,kk}(:,3)).^2 + sdata{jj,kk}(:,1).^2./xdata{jj,kk}(:,1).^2);
        sig_ihMTR = sqrt(var_ihMTR);
        %%% 95% confidence interval - valid to do here rather than on
        %%% individual variables because n is the same for all
        ci_ihMTR = tinv([0.975],sdata{jj,3}-1) * sig_ihMTR / sqrt(sdata{jj,3});
        
        MTR = 100*(xdata{jj,kk}(:,1)-xdata{jj,kk}(:,2))./xdata{jj,kk}(:,1);
        var_MTR = MTR.^2 .* ( (sdata{jj,kk}(:,1).^2+sdata{jj,kk}(:,2).^2)./(xdata{jj,kk}(:,1)-xdata{jj,kk}(:,2)).^2 + sdata{jj,kk}(:,1).^2./xdata{jj,kk}(:,1).^2);
        sig_MTR = sqrt(var_MTR);
        %%% 95% confidence interval
        ci_MTR = tinv([0.975],sdata{jj,3}-1) * sig_MTR / sqrt(sdata{jj,3});
        
        
        %%% Plot Data
        if kk==1
            tmp3=errorbar(seq_pars.flips*180/pi,MTR,ci_MTR,'o','markersize',ms,'color',[0 0 1]*0.5,'markerfacecolor',[0.5 0.5 1]);
            tmp1=errorbar(seq_pars.flips*180/pi,ihMTR,ci_ihMTR,'o','markersize',ms,'color',[0 1 0]*0.5,'markerfacecolor',[0.5 1 0.5]);
        else
            tmp4=errorbar(seq_pars_spgr.flips*180/pi,MTR,ci_MTR,'^','markersize',ms,'color',[0 0 1]*0.5,'markerfacecolor',[0.5 0.5 1]);
            tmp2=errorbar(seq_pars_spgr.flips*180/pi,ihMTR,ci_ihMTR,'^','markersize',ms,'color',[0 1 0]*0.5,'markerfacecolor',[0.5 1 0.5]);
        end
    end
    grid on
    grid on
    xlabel('Flip angle, degrees')
    ylabel('Ratio (% of S_{1B})')
    
    if jj==1
        ll=legend('MTR, bSSFP','ihMTR, bSSFP','MTR, SPGR','ihMTR, SPGR','AutoUpdate','off');
        ll=legend('MTR, bSSFP','MTR, SPGR','ihMTR, bSSFP','ihMTR, SPGR','AutoUpdate','off');
        %ll.Location = 'EastOutside';
        ll.FontSize = 12;
    end
    yl = get(gca,'YLim');
    ylim([yl(1) 60])
    
    gg=gca;
    gg.XAxisLocation = 'origin';

    title(titles{jj},'fontsize',12)
end
setpospap([100 100 850 550])


%% Repeat above but with proper error propagation for MTR and ihMTR
%% Second version, with separate y-axes

titles = {'Water with MnCl_2','BSA','PL-161','Hair Conditioner'};
nr=2;nc=2;
ms=4;

figfp(1)

for jj=1:4
   
    subplot(nr,nc,jj)
    
    hold on 
    grid on
    
    for kk=1:2 %<--- loop over sequence type 1=bSSFP, 2=SPGR
        
        ihMTR = 100*(xdata{jj,kk}(:,2)-xdata{jj,kk}(:,3))./xdata{jj,kk}(:,1);
        var_ihMTR = ihMTR.^2 .* ( (sdata{jj,kk}(:,2).^2+sdata{jj,kk}(:,3).^2)./(xdata{jj,kk}(:,2)-xdata{jj,kk}(:,3)).^2 + sdata{jj,kk}(:,1).^2./xdata{jj,kk}(:,1).^2);
        sig_ihMTR = sqrt(var_ihMTR);
        %%% 95% confidence interval - valid to do here rather than on
        %%% individual variables because n is the same for all
        ci_ihMTR = tinv([0.975],sdata{jj,3}-1) * sig_ihMTR / sqrt(sdata{jj,3});
        
        MTR = 100*(xdata{jj,kk}(:,1)-xdata{jj,kk}(:,2))./xdata{jj,kk}(:,1);
        var_MTR = MTR.^2 .* ( (sdata{jj,kk}(:,1).^2+sdata{jj,kk}(:,2).^2)./(xdata{jj,kk}(:,1)-xdata{jj,kk}(:,2)).^2 + sdata{jj,kk}(:,1).^2./xdata{jj,kk}(:,1).^2);
        sig_MTR = sqrt(var_MTR);
        %%% 95% confidence interval
        ci_MTR = tinv([0.975],sdata{jj,3}-1) * sig_MTR / sqrt(sdata{jj,3});
        
        
        %%% Plot Data
        if kk==1
            yyaxis left
            tmp3=errorbar(seq_pars.flips*180/pi,MTR,ci_MTR,'o','markersize',ms,'color',[0 0 1]*0.5,'markerfacecolor',[0.5 0.5 1]);
            yyaxis right
            tmp1=errorbar(seq_pars.flips*180/pi,ihMTR,ci_ihMTR,'o','markersize',ms,'color',[0 1 0]*0.5,'markerfacecolor',[0.5 1 0.5]);
        else
            yyaxis left
            tmp4=errorbar(seq_pars_spgr.flips*180/pi,MTR,ci_MTR,'^','markersize',ms,'color',[0 0 1]*0.5,'markerfacecolor',[0.5 0.5 1]);
            yyaxis right
            tmp2=errorbar(seq_pars_spgr.flips*180/pi,ihMTR,ci_ihMTR,'^','markersize',ms,'color',[0 1 0]*0.5,'markerfacecolor',[0.5 1 0.5]);
        end
    end
    grid on
    grid on
    xlabel('Flip angle, degrees')
    
    
    if jj==1
        ll=legend('MTR, bSSFP','MTR, SPGR','ihMTR, bSSFP','ihMTR, SPGR','AutoUpdate','off');
        %ll.Location = 'EastOutside';
        ll.FontSize = 12;
    end
    yyaxis left
    yl = get(gca,'YLim');
    ylim([-5 60])
    ylabel('MTR (% of S_{1B})')
    gg=gca;
    gg.YColor = [0 0 1]*0.5;
    yyaxis right
    yl = get(gca,'YLim');
    ylim([-20/12 20])
    ylabel('ihMTR (% of S_{1B})')
    gg=gca;
    gg.YColor = [0 1 0]*0.5;
    %gg.XAxisLocation = 'origin';

    title(titles{jj},'fontsize',12)
end
setpospap([100 100 850 550])


%% Bootstrapping for uncertainty determination

% store results here
fitpars = {};
fitdata = {};
fiterr =  {};

%%% Number of trials
Ntrial = 500;


%%% SELECT VOXEL AND PHANTOM %%%%%%%
for ix = 1:4 %<--- Loop over phantom
    

    fitpars{ix} = zeros(Ntrial,10);
    fitdata{ix} = zeros(8,6,Ntrial);
    fiterr{ix} = zeros(Ntrial,1);
    
    % Get the GROUND TRUTH solution and forward model
    x0 = sol{ix};
    s0 = fwd_model(x0);
    
    data0 = cat(2,xdata{ix,1},xdata{ix,2});
    residual0 = fwd_model(x0)-data0;
    
    % For each trial, augment data by adding permuted residuals
    for ii=1:Ntrial
 
        pix = ceil(48*rand(48,1)); %<-- we want to randomize with replacement, not permute
        
        residual_perm = residual0;
        residual_perm(:) = residual_perm(pix);
        
        data_mod = s0 + residual_perm;
        
        costf = @(x)(norm(fwd_model(x)-data_mod)^2);
               
        options = optimoptions('fmincon');
        options.FiniteDifferenceStepSize = 1e-4; %<-- use coarser finite difference steps to speed convergence
              
        [xs,xf,xflg,output] = fmincon(costf,x0,[],[],[],[],lb{ix},ub{ix},[],options);
        
        fitpars{ix}(ii,:)=xs;
        
        %%% Store forward prediction
        fitdata{ix}(:,:,ii) = fwd_model(xs);
    
        %%% Store fit error wrt original data
        fiterr{ix}(ii) = 100*norm(data0-fitdata{ix}(:,:,ii))/norm(data0);
        
        fprintf(1,'Phantom %d of 4, trial number %d of %d\n',ix,ii,Ntrial);
    end
end

save bin/bootstrapdata fitpars fitdata fiterr


%% Print data


names = {'MnCl2','BSA','PL161','HC'};

for ix=1:4
    
    
    %%% Remove outliers for HC fitting
    if ix==4
        idx = fiterr{ix}<5;
        xm = mean(fitpars{ix}(idx,:),1);
        xs = std(fitpars{ix}(idx,:),1);
    else
        xm = mean(fitpars{ix},1);
        xs = std(fitpars{ix},1);
    end
    fprintf(1,'\n%s **** \n  R1f  = %1.3f ± %1.3f s^-1 \n  R2f  = %1.3f ± %1.3f s^-1 \n  R1SZ = %1.3f ± %1.3f s^-1 \n  R1D  = %1.3f ± %1.3f s^-1 \n  K  = %1.3f ± %1.3f s^-1\n',...
        names{ix},xm(1),xs(1),xm(2),xs(2),xm(4),xs(4),xm(5),xs(5),xm(8),xs(8))
    
    fprintf(1,'\n%s **** \n  T1f  = %1.0f ± %1.0f ms \n  T2f  = %1.1f ± %1.1f ms \n  T1SZ = %1.1f ± %1.1f ms \n  T1D  = %1.1f ± %1.1f ms \n  T2s  = %1.1f ± %1.1f us \n',...
        names{ix},1000/xm(1),1000*xs(1)/xm(1)^2,1000/xm(2),1000*xs(2)/xm(2)^2,1000/xm(4),1000*xs(4)/xm(4)^2,1000/xm(5),1000*xs(5)/xm(5)^2,xm(7),xs(7))
    
    fprintf(1,'  M0s  = %1.3f ± %1.3f \n  f    = %1.3f ± %1.3f \n',xm(3),xs(3),xm(6),xs(6))
    
end


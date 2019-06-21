%%% Fitting script bootstrap error estimate

addpath(genpath('../../../published code/ihMT_steadystate_MB'))

%% Define SSFP pulse sequence
seq_pars = struct;
seq_pars.flips = d2r(10:10:80);
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

% Define SPGR pulse sequence
seq_pars_spgr = seq_pars;
seq_pars_spgr.flips = d2r(2:2:16);
nf = length(seq_pars_spgr.flips);

seq_pars_spgr.b1sqrd = {};

for ii=1:nf
    % 3 bands
    [~,seq_pars_spgr.b1sqrd{ii,3},~] = gen_MB_pulse(seq_pars_spgr.flips(ii),seq_pars_spgr.tau,seq_pars_spgr.TR,seq_pars_spgr.b1_rms,seq_pars_spgr.delta,'3','sigma',3,'dt',dt); 
    % 2 bands
    [~,seq_pars_spgr.b1sqrd{ii,2},~] = gen_MB_pulse(seq_pars_spgr.flips(ii),seq_pars_spgr.tau,seq_pars_spgr.TR,seq_pars_spgr.b1_rms,seq_pars_spgr.delta,'2+','sigma',3,'dt',dt);
    % 1 bands
    seq_pars_spgr.b1sqrd{ii,1} = seq_pars_spgr.b1sqrd{ii,2}.*[0 1 0]; 
end
seq_pars_spgr.Delta_Hz = seq_pars.delta * [-1 0 1];

% add echo time
seq_pars_spgr.TE = 2.5e-3;

%% Load in data to fit
load fitdata.mat
IX = 3; %PL161
IX = 2;% BSA
IX = 4; % Hair cond
IX = 5; % MnCl2

IX = [5 2 3 4];

%%% Get the data for each phantom into a cell array
xdata = {};

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
    
    % scale down a bit
    xdata{jj,1} = permute(tmp/1000,[2 3 1]);
    xdata{jj,2} = permute(tmp2/1000,[2 3 1]);
    
    %%% For bootstrapping look at means
    xdata{jj,1} = mean(xdata{jj,1},3);
    xdata{jj,2} = mean(xdata{jj,2},3);
end


%% combined cost function and fit

%%% Solutions from mean fit - these are starting points for voxel-wise
% xs_HC = [0.6039    6.4390    0.0613    3.1649   43.6045    0.7880   20.9714   77.8662    5.0962    1.0000]; %<--- not fixed
% xs_BSA= [0.5849   19.4224    0.0849    6.7534   46.8443         0   12.6667   62.7665    4.8752    1.0000];
% xs_PL = [0.5513   14.0496    0.1517    4.4166   48.2816    1.0000   17.8180   46.6658    4.9998    1.0000];
% xs_Mn = [0.5206    1.6642         0    1.2068   49.9994         0   12.0212   50.0003    5.2583    1.0000];
% Updated solutions, here because TBP was altered in Bieri-Scheffler
% correction
xs_HC = [0.5955    5.9999    0.0675    2.9780   42.8042    0.7139   25.8479   67.8720    5.1397    1.0000];
xs_BSA= [0.5635   18.8122    0.0774    7.7545   47.8069         0   15.2003   70.6874    4.8311    1.0000];
xs_PL = [0.5515   13.6859    0.1516    4.4111   48.2443    1.0000   17.8722   46.7160    4.9972    1.0000];%<-- fix f=1
xs_Mn = [0.5205    1.6267         0    4.8914   48.2782    0.4920   17.8562   46.6711    5.2583    1.0000];

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

    
sol = {};
sol{1} = xs_Mn;
sol{2} = xs_BSA;
sol{3} = xs_PL;
sol{4} = xs_HC;

%%% Store results
fitpars = {};
fitdata = {};
fiterr =  {};

%%% Number of trials
Ntrial = 500;


%%% Single forward model function
flag = 1; %<-- Gaussian line
fwd_model = @(x)(cat(2,ssfp_ihmt_fit_fwd_model(x,seq_pars,flag),spgr_ihmt_fit_fwd_model(x,seq_pars_spgr,flag)));

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
 
        %pix = randperm(48);
        pix = ceil(48*rand(48,1)); %<-- we want to randomize with replacement, not permute
        
        residual_perm = residual0;
        residual_perm(:) = residual_perm(pix);
        
        data_mod = s0 + residual_perm;
        
        costf = @(x)(norm(fwd_model(x)-data_mod)^2);
               
        options = optimoptions('fmincon');
        %options = optimoptions(options,'Display', 'iter');
        options.FiniteDifferenceStepSize = 1e-4; %<-- use coarser finite difference steps to speed convergence
        %options = optimoptions(options,'PlotFcns', {  @optimplotx @optimplotfunccount @optimplotfval });
        
        [xs,xf,xflg,output] = fmincon(costf,x0,[],[],[],[],lb{ix},ub{ix},[],options);
        
        fitpars{ix}(ii,:)=xs;
        
        %%% Store forward prediction
        fitdata{ix}(:,:,ii) = fwd_model(xs);
    
        %%% Store fit error wrt original data
        fiterr{ix}(ii) = 100*norm(data0-fitdata{ix}(:,:,ii))/norm(data0);
        
        fprintf(1,'Phantom %d of 4, trial number %d of %d\n',ix,ii,Ntrial);
    end
end

%save bootstrapdata2 fitpars fitdata fiterr
%save bootstrapdata3 fitpars fitdata fiterr
% save bootstrapdata4 fitpars fitdata fiterr


%% Print data


names = {'MnCl2','BSA','PL161','HC'};

for ix=1:4
    
    
    %%% Remove outliers for HC fitting
    if ix==4
        idx = fiterr{ix}<4;
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

%% 
for ix=1:4
    figfp(ix)
    for ii=1:10
        subplot(3,4,ii)
        histogram(fitpars{ix}(:,ii))
    end
end
figfp(5)

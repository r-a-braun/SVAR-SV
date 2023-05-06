clear; clc; rng(100)
addpath('functions_svsvar')
%% This code illustrates the use of the model and the computational trade-off between EM-1 and EM-2

%% Simulate a SV-SVAR model
A = [0,0.6,0.35; 0,-0.1,.7]; %intercept and VAR slopes
B = [1 0;.5 2]; % impact matrix
phi = .95; % persistance parameters
s = .04; % shock variance parameters
mink = 3.6; % minimum kurtosis
p = 1; % lag length
T = 200; % time series length
bi = 50; % burn-in
rep = 1; % no. of replications
y_sim = Simulate_SV(A,B,phi,s,T,bi,rep,mink);
y = y_sim(:,:,1)';

%% Estimate Models and compare 
[T, K] = size( y );
specs.p = 1;  % lag length
specs.max_it = 5000;  % Maximum # iterations EM algorithm
specs.M = 20000; % Draws Importance Sampler Likelihood evaluation
specs.M2 = 20000; % Draws Importance Sampler EM algorithm
specs.tol = 1e-04; % Tolerance EM algorithm
specs.alpha = .1; % Confidence level IRFs/FEVDs
specs.hor = 16; % Horizons IRF/FEVDs
specs.S_b = eye(K^2); % selection matrix B elements (unrestricted)
specs.scale = 20; % offset constant : std(error)./scale
specs.r = K; % number of heteroskedastic shocks
specs.unitshocks  = 1; % if 1, computes IRFs to shocks of size 1, else s.d.
specs.param_init = [];

% EM-1 parameter estimates only (~3 sec)
tic;
[Param_est_EM1 ] = SVAR_SV_em_approx( y, specs);
toc;
disp([ B, Param_est_EM1.B./diag(Param_est_EM1.B)'])

% EM-1 parameter plus likelihood, covariance matrix and IRFs/FEVDs (~ 5 sec)
tic;
[Param_est_EM1, loglike_EM1, count_EM1, VarTheta_EM1, IRF_EM1] = SVAR_SV_em_approx( y, specs);
toc; 

% EM-2 parameter plus likelihood, covariance matrix and IRFs/FEVDs (~ 1 minute)
tic;
[Param_est_EM2, loglike_EM2, count_EM2, VarTheta_EM2, IRF_EM2]  = SVAR_SV_em_IS( y, specs);
toc; 

%% Compare IRFs: EM1 (blue) vs EM2 (red)
IRF_true = IRFvar( A, B , specs.hor, 1); 
variables = {' y1',  ' y2'};
shocks = {' e1',  ' e2'};
horizon = (0:specs.hor)';
a = 1;
for i = 1:K
    for ii = 1:K % last shock corresponds to MP
        subplot(K,K,a)  
        plot(horizon, IRF_true.irfs(:,(ii*K-K)+i ),'k',...
            horizon, IRF_EM1.irfs(:,(ii*K-K)+i ),'b.-', ...
            horizon,[squeeze(IRF_EM1.IRF_qu( :, (ii*K-K)+i,:)), squeeze(IRF_EM1.IRF_ql(:,(ii*K-K)+i,:))],'b--',...
            horizon, IRF_EM2.irfs(:,(ii*K-K)+i ),'r.-', ...
            horizon,[squeeze(IRF_EM2.IRF_qu( :, (ii*K-K)+i,:)), squeeze(IRF_EM2.IRF_ql(:,(ii*K-K)+i,:))],'r--')
        xlim([horizon(1),horizon(end)])
        title(strcat('$',shocks{i},'\rightarrow',variables{i},'$'),'Interpreter','latex')
        a = a + 1;
    end
end


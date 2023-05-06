clear; clc;
load('data/BH19.mat')
addpath('functions_svsvar')
y = [100*diff(log(dataBH19.production)),...
    100*diff(log(dataBH19.OECD)),...
    100*log(dataBH19.realWTI(2:end)),...
    dataBH19.dInvent(2:end)];
date = dataBH19.date(2:end); 
startdate = datetime(1974,1,1);  
idx_Y = datenum(date)>=datenum(startdate); 
y = y(idx_Y,:); 
[T,K]=size(y);

 %% Estimate Model via EM-2 (Importance Sampling) 
specs.p = 12;  % lag length
specs.max_it = 5000;  % Maximum # iterations EM algorithm
specs.M = 20000; % Draws Importance Sampler Likelihood evaluation
specs.M2 = 20000; % Draws Importance Sampler EM algorithm  
specs.tol = 1e-05; % Tolerance EM algorithm
specs.alpha = [.1, .32]; % Confidence level IRFs/FEVDs
specs.hor = 72; % Horizons IRF/FEVDs
specs.S_b = eye(K^2); % selection matrix B elements (unrestricted)
specs.scale = 20; % offset constant : std(error)./scale
specs.r = K; % number of heteroskedastic shocks
specs.unitshocks  = 1; % if 1, computes IRFs to shocks of size 1, else s.d.
tic;
[Param_est , loglike, count, VarTheta, IRF]  = SVAR_SV_em_IS( y, specs);
toc; 
save('output/results_bh19.mat')
 
% %% Plot IRFs SV-SVAR model
variables = {' prod',  ' Activity', ' price', ' invent'};
shocks = {' e1',  ' e2', ' e3', ' e4'};
csdummy = [1,1,0,1];
horizon = (0:specs.hor)';  
a = 1;
for i = 1:K
    for ii = 1:K % last shock corresponds to MP
        subplot(K,K,a)
        if csdummy(i)==1
            plot(horizon, IRF.irfs_cs(:,(ii*K-K)+i ),'k.-', ...
            horizon,[squeeze(IRF.IRF_cs_qu( :, (ii*K-K)+i,:)), squeeze(IRF.IRF_cs_ql(:,(ii*K-K)+i,:))],'k.-')
        else
            plot(horizon, IRF.irfs(:,(ii*K-K)+i ),'k.-', ...
            horizon,[squeeze(IRF.IRF_qu( :, (ii*K-K)+i,:)), squeeze(IRF.IRF_ql(:,(ii*K-K)+i,:))],'k.-') 
        end
        
        xlim([horizon(1),horizon(end)])
        title(strcat('$',shocks{i},'\rightarrow',variables{i},'$'),'Interpreter','latex')
        a = a + 1;
    end
end
% 
% 

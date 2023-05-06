clear; clc; 
%% This code creates Table 1
% Replication code for OTHER models is not included in the file
load('data/BH19.mat')
y = [100*diff(log(dataBH19.production)),...
    100*diff(log(dataBH19.OECD)),...
    100*log(dataBH19.realWTI(2:end)),...
    dataBH19.dInvent(2:end)];
date = dataBH19.date(2:end); 
startdate = datetime(1974,1,1);  
idx_Y = datenum(date)>=datenum(startdate); 
y = y(idx_Y,:);  
[T, K]=size(y); 
p = 12;
[Bhat, Sigmahat, Uhat, llh_lin] = VARmlike(y,p) ;
n_param = K*(K*(p+1)) + K*(K+1)/2; 
AIC_lin = -2*llh_lin + 2*n_param;
BIC_lin = -2*llh_lin + log(T)*n_param;
Data(:,1) = [llh_lin;AIC_lin;BIC_lin];
load('output/Estimate_other_models/MS2.mat')  
Data(:,2) = [LLH_MS2.llh,LLH_MS2.AIC,LLH_MS2.BIC];
load('output/Estimate_other_models/MS3.mat')  
Data(:,3) = [LLH_MS3.llh,LLH_MS3.AIC,LLH_MS3.BIC];
load('output/Estimate_other_models/STVAR.mat')  
Data(:,4) = [LLH_STVAR.llh,LLH_STVAR.AIC,LLH_STVAR.BIC];
load('output/Estimate_other_models/GARCH.mat')  
Data(:,5) = [LLH_GARCH.llh,LLH_GARCH.AIC,LLH_GARCH.BIC]; 
load('output/results_bh19.mat') 
Data(:,6) = [loglike.llh,loglike.aic,loglike.bic];
input.data = Data;
input.tableColLabels = {'linear', 'MS($2$)', 'MS($3$)', 'STVAR', 'GARCH($1,1$)', 'SV(1)'};
input.tableRowLabels = {'log likelihood','AIC','BIC'};
input.dataFormat = {'%.2f',6}; % three digits precision for first two columns, one digit for the last
input.tableCaption = 'Model selection';
input.tableLabel = 'modelselec';
latex = latexTable(input); 

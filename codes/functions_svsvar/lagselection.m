function [ FPE, AIC, HQ, BIC] = lagselection(Yraw,pmax, inc) %Function gives the lag order chosen by the criterias

FPEraw = ones(pmax+1,1);
AICraw = ones(pmax+1,1);
HQraw = ones(pmax+1,1);
BICraw = ones(pmax+1,1);
[Tpmax, K] = size(Yraw);
T = Tpmax - pmax;
for m = 0:pmax 
    Y = Yraw((1+(pmax-m)):end,:); % same estimation length for all
    [~,Sigmahat] = VARmvls(Y,m,inc); 
    Sigmahat = ((T-K*m-inc)/T)*Sigmahat; % Max Like Variance Estimator
    FPEraw(m+1) = ((T+K*m+1)/(T-K*m-1))^K * det(Sigmahat); %Criteria Values
    AICraw(m+1) = log(det(Sigmahat)) + (2*m*K^2)/T;
    HQraw(m+1) = log(det(Sigmahat)) + (2* log(log(T)))/T * m * K^2;
    BICraw(m+1) = log(det(Sigmahat)) + (log(T))/T * m* K^2;
end
[~,FPE]=min(FPEraw);  
FPE=FPE-1;
[~,AIC]=min(AICraw);
AIC=AIC-1;
[~,HQ]=min(HQraw);
HQ=HQ-1;
[~,BIC]=min(BICraw);
BIC=BIC-1;

end 

function [Bhat, Sigmahat, Uhat, AValpha, AVsigma] = VARmvls(y,p,inc)
%VAR Computes Maximum Likelihood estimates of VAR(p) with intercept (not
%mean adjusted)
% inputs:   -y: the data of size (TxK)
%           -p: VAR lag order,
%           -inc: intercept dummy
% outputs:  -Bhat: ML estimates AR parameters
%           -Sigmahat: ML estimate COV matrix
%           -Uhat: residuals
%           -Tstatistic: t stats of parameters

[TpP, K] = size(y);
T = TpP - p;
Y = y(p+1:end,:)';  
rhs = lagmatrix(y,1:p);
Z = rhs(p+1:end,:)';
if inc == 1
    Z = [ones(1,T); Z];
end
Bhat = (Y*Z')/(Z*Z'); %LS Estimator
Uhat = (Y-Bhat*Z);
Sigmahat = 1/(T-K*p-inc)*(Uhat)*Uhat'; % Residual Cov-Matrix 
%Sigmahat = 1/(T)*(Uhat)*Uhat'; % Residual Cov-Matrix  
if nargout>3 
    AValpha = kron(inv(Z*Z'/T),Sigmahat); % Cov Matrix of Estimated Parameters
    D = duplication(K);
    Dk = (D'*D)\D';
    AVsigma = 2*Dk*kron(Sigmahat,Sigmahat)*Dk';
end
end
 

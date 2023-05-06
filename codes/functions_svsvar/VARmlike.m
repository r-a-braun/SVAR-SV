function [Bhat, Sigmahat, Uhat,llh] = VARmlike(y,p)
%VAR Computes LS estimates of VAR(p) with intercept
% inputs:   -y: the data of size (TxK)
%           -p: VAR lag order,
%           
% outputs:  -Bhat: OLS estimates AR parameters
%           -Sigmahat: OLS estimate COV matrix
%           -Uhat: residuals
%          
[TpP, K] = size(y);
T = TpP - p;
Y = y(p+1:end,:)';  
rhs = lagmatrix(y,1:p);
Z = rhs(p+1:end,:)';
Z = [ones(1,T); Z]; %add the one's on top, to include the intercept

Bhat = Y*Z'*inv(Z*Z'); %LS Estimator
Uhat = (Y-Bhat*Z)';
Sigmahat = 1/T*Uhat'*Uhat; % Residual Cov-Matrix

if nargout>3 
     llh = -T/2*K*log(2*pi) - T/2*log(det(Sigmahat)) - 1/2*trace(Uhat*inv(Sigmahat)*Uhat');
end


 
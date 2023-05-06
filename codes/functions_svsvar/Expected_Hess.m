function [H] = Expected_Hess(param_est,Y,X,H_i,Var_H,Cov_Htm1,iV,S_b)
% Bertsche, Braun (2020): Identification of Structural Vector
% Autoregressions by Stochastic Volatility
%
% function to compute gradients and hessians required for Louis identity
%
%% Input:
% - para (2r + m*K + K^2 x 1): [phi;s;vec(A);vec(B)]
% - Y (T x K): dependent variables
% - X (T x Kp+ic): explanatory variables
% - Eh (T x K x nrep): expected log variances
% - Covh (T x T x r x nrep): covariance matrices of log variances
% - S_B (K^2 x nB): selection matrix for vec(B)
%
%% Output:
% - g2out (npara x npara x nrep): outer product of gradients
% - g2in (npara x npara x nrep): inner product of gradients
% - H (npara x npara x nrep): hessian

%% Preparations
[T,K] = size(Y); %no. of time points, heteroskedastic shocks and replications
m = size(X,2); %no. of parameters per VAR equation 
phi = param_est.phi;  %AR parameters SV equation
s = param_est.s;  %variance parameters SV equation 
A =  param_est.A;   %VAR parameters 
B =  param_est.B;  %impact matrix
iB = inv(B); %inverse impact matrix 
K_KK = commutation(K,K); %commutation matrix
Ik = eye(K); %identity
KpIk2 = K_KK+eye(K.^2); %expression required later 
u = Y-X*A; %reduced form errors
eps = (B\u')'; %structural shocks  
idx_r = find(isnan(phi)==0);
r = length(idx_r);
npara = length([phi(idx_r);s(idx_r);vec(A);pinv(S_b)*vec(B)]);

%% Prepare empty matrices 
H = zeros(npara,npara); 
vSxx = sum(Var_H(1:end-1,:) + H_i(1:end-1,:).^2)';
vSxy = sum( Cov_Htm1+H_i(1:end-1,:).*H_i(2:end,:) )';
vSyy = sum(Var_H(2:end,:)+H_i(2:end,:).^2)';
Eh1sq = (Var_H(1,:) + H_i(1,:).^2)'; 
Hess_corr = zeros(2*r);
for i = idx_r'
    [~, Hess_ci] = gradhess_correction( [phi(i);s(i)] , T); 
    Hess_corr(i:r:end,i:r:end) = Hess_ci; 
end    
Hphi = (- 1./s.*vSxx ) + (-(1+phi.^2)./((1-phi.^2).^2)  + Eh1sq./s );
Hphis= (-1./(s.^2).*vSxy + 1./(s.^2).*vSxx.*phi) +  (-phi./(s.^2).*Eh1sq);
Hs = ( (T-1)./2*(1./(s.^2)) - 1./(s.^3).*(vSyy - 2*vSxy.*phi  + phi.^2.*vSxx) )...
    + (1./(2.*s.^2)-(1-phi.^2)./(s.^3).*Eh1sq); 
H(1:2*r,1:2*r) = diag([Hphi(idx_r);Hs(idx_r)])+diag(Hphis(idx_r),-r)+diag(Hphis(idx_r),r) + Hess_corr;
 
%% Gradient and hessian in A and B 
ScoreBeta = zeros(K); %to be filled for gradient of beta
Hess_alpha = zeros(m*K); %empty matrix for hessian wrt alpha
Hess_AB = zeros(m*K,K^2); %empty matrix for cross derivatives
Hess_beta = T*kron(iB,iB')*K_KK; %empty matrix for hessian wrt beta
for t = 1:T
    iVt = diag(iV(t,:)); %inverse variance of structural shocks
    Xtil = kron(X(t,:),Ik); %matrix form of explanatory variables
    iSigt = iB'*iVt*iB; %inverse variance of reduced form errors
    iSXtil = iSigt*Xtil; %matrix product required several times
    SumScore = Xtil'*iSXtil; %matrix product required several times
    Hess_alpha = Hess_alpha-SumScore; %sum of 2nd derivative wrt alpha
    Eps2 = eps(t,:)'*eps(t,:); %squared structural shocks 
    ScoreBeta = ScoreBeta+iB'*iVt*Eps2; %sum for gradient wrt beta
    Hess_AB = Hess_AB-(kron(eps(t,:),iSXtil')+kron(u(t,:)*iSigt,Xtil'*iB')*K_KK); %cross derivative wrt alpha and beta
    Hess_beta = Hess_beta-kron(Ik,iB'*iVt)*KpIk2*kron(Eps2,iB)-...
        kron(Eps2*iVt*iB,iB')*K_KK; %2nd derivative wrt to beta
end 
H(2*r+1:2*r+m*K,2*r+1:2*r+m*K) = Hess_alpha;
H(2*r+m*K+1:end,2*r+1:2*r+m*K) = (Hess_AB*S_b)';
H(2*r+1:2*r+m*K,2*r+m*K+1:end) = Hess_AB*S_b;
H(2*r+m*K+1:end,2*r+m*K+1:end) = S_b'*Hess_beta*S_b;

end

function a = vec(A)
a = A(:);
end
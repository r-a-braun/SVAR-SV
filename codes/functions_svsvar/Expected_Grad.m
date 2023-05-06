function g = Expected_Grad(para_exp,Param_est,Y,X,p, h_lap,S_B,scale,CRNo)
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
m = size(X,2); %no. of parameters per VAR equat  
phi = Param_est.phi;
s = Param_est.s; 
A = Param_est.A; %VAR parameters 
B = Param_est.B; %impact matrix
iB = inv(B); %inverse impact matrix
va = vec(A'); %vectorized transposed A  
Ik = eye(K); %identity  
u = Y-X*A; %reduced form errors
eps = (B\u')'; %structural shocks
npara = length(para_exp); %no. of parameters 
idx_r = find(isnan(phi)==0)';  
r = length(idx_r); 
% Compute Expectations: 
nA = K*(K*p+1);  
phie = NaN(K,1); se = NaN(K,1);
phie(idx_r) = para_exp(1:r); 
se(idx_r) = para_exp(r+1:2*r);     
vAlpe = para_exp(2*r+1:2*r+nA);
vBe = para_exp(2*r+nA+1:end);  
Ae = reshape(vAlpe,K,(K*p+1))'; 
Be = reshape(S_B*vBe,K,K); 
ee = (Be\(Y-X*Ae)')';  

%% E-STEP
H_i = zeros(T,K); Var_H = zeros(T,K); Cov_Htm1 = zeros(T-1,K); iV= ones(T,K);
for j = idx_r 
    if nargin>8
        [H_i(:,j), Var_H(:,j), Cov_Htm1(:,j),iV(:,j) ] = e_step_IS(ee(:,j), phie(j), se(j), h_lap(:,j), CRNo(:,:,j) , scale  );
    else
        [H_i(:,j), Var_H(:,j), Cov_Htm1(:,j),iV(:,j)] = e_step_approx(ee(:,j), phie(j), se(j), h_lap(:,j) , scale);
    end  
end
    
     
%% Gradient phi s
grad_corr = zeros(2*r,1); 
for i = idx_r
    grad_ci = gradhess_correction( [phi(i);s(i)] , T);
    grad_corr(i:r:end) = grad_ci;  
end 
g = zeros(npara,1); %empty matrix for gradient
vSxx = sum(Var_H(1:end-1,:) + H_i(1:end-1,:).^2)';
vSxy = sum( Cov_Htm1+H_i(1:end-1,:).*H_i(2:end,:) )';
vSyy = sum(Var_H(2:end,:)+H_i(2:end,:).^2)';  
Eh1sq = (Var_H(1,:) + H_i(1,:).^2)';
gphi1 = -phi./(1-phi.^2) + phi./s.*Eh1sq;
gphi = 1./(s).*vSxy - 1./s.*vSxx.*phi  + gphi1;
gs1 = - 1./(2.*s) + ((1-phi.^2)./(2.*s.^2)).*Eh1sq;
gs = -(T-1)./2*(1./(s)) + 1./(2.*s.^2).*(vSyy - 2*vSxy.*phi  + phi.^2.*vSxx)  + gs1; 
g(1:2*r) = [gphi(idx_r);gs(idx_r)]+grad_corr;  


%% Gradient in A and B 
grad_alpha = zeros(m*K,1); %empty matrix for gradient wrt alpha
ScoreBeta = zeros(K); %to be filled for gradient of beta 
for t = 1:T
    iVt = diag(iV(t,:)); %inverse variance of structural shocks
    Xtil = kron(X(t,:),Ik); %matrix form of explanatory variables
    iSigt = iB'*iVt*iB; %inverse variance of reduced form errors
    iSXtil = iSigt*Xtil; %matrix product required several times
    SumScore = Xtil'*iSXtil; %matrix product required several times 
    Eps2 = eps(t,:)'*eps(t,:); %squared structural shocks
    grad_alpha = grad_alpha+((Y(t,:)*iSigt*Xtil)-va'*SumScore)'; %1st derivative wrt alpha
    ScoreBeta = ScoreBeta+iB'*iVt*Eps2; %sum for gradient wrt beta 
end
g(2*r+1:2*r+m*K) = grad_alpha; %gradient wrt A
g(2*r+m*K+1:end) = ((-T*vec(iB')'+vec(ScoreBeta)')*S_B)'; %gradient wrt B

end

function [ Q1,Q2,LM ] = svar_sv_tests( uhat, Bhat, r, Hw)
%SVAR_TESTS_IDENTIFICATION Summary of this function goes here
%   Detailed explanation goes here
A = inv(Bhat); 
A2 = A(r+1:end,:);
[T,n] = size(uhat); 
xi = zeros(T,1);
n_nu = (n-r)*(n-r+1)/2;
eta = zeros(T,n_nu);

for t=1:T
    xi(t) = uhat(t,:)*(A2)'*A2*uhat(t,:)';
    eta(t,:) = vech(A2*uhat(t,:)'*uhat(t,:)*A2'); 
end
xi = xi-mean(xi);
nu = eta - mean(eta);

% Q1 Test of Lanne and Saikkonen 
gamma = zeros(Hw,1);
gamma0 = 1/T*sum(xi.^2);
for h = 1:Hw
    gamma(h) = 1/T*sum(xi(1:end-h).*xi(h+1:end));
end 
Q1.teststat = T*sum((gamma./gamma0).^2);
Q1.dof = Hw;
Q1.p = 1 - chi2cdf(Q1.teststat,Q1.dof); 

% Q2 Test of Lanne and Saikkonen 
G0 = 1/T*(nu')*nu;
iG0 = inv(G0);
Q2s = 0;
for h = 1:Hw
    Gh = 1/T*(nu(1:end-h,:)'*nu(h+1:end,:));
    Q2s = Q2s + T*trace(Gh'*iG0*Gh*iG0);
end
Q2.teststat = Q2s;
Q2.dof = 1/4*Hw*(n-r)^2*(n-r+1)^2;
Q2.p = 1 - chi2cdf(Q2.teststat,Q2.dof); 

% Lütkepohl and Milunovich Test based on LM for residuals
[~, Sig_ur] = VARmlike(eta,Hw);
LM.teststat = 1/2*T*(n-r)*(n-r+1)-T*trace(Sig_ur*iG0);
LM.dof = 1/4*Hw*(n-r)^2*(n-r+1)^2;
LM.p = 1 - chi2cdf(LM.teststat,LM.dof); 

 

end


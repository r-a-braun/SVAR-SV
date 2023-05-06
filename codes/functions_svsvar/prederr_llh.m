function [ llh ] = prederr_llh(Y, X ,p, phi, s, vApvB )
% Approximate log likelihood of the SV-SVAR model based on prediction error decomposition
 
%% Read Specifications:  
[T,K] = size(Y);    
nA = K*(K*p+1);   
vAlp = vApvB(1:nA);
vB = vApvB(nA+1:end);  
A = reshape(vAlp,K,(K*p+1))'; 
B = reshape(vB,K,K);  
U = Y - X*A;
Err = (B\U')';  
s2 = Err.^2; 
 
%% Filter algorithm to evaluate log likelihood for each shock
llh_t = zeros(T,K);
hfil = zeros(T,K); 
%gamma2s = zeros(T,K); betats = zeros(T,K); sigma2s = zeros(T,K);
for i = 1:K 
    betat = 0;
    gamma2 = s(i)/(1-phi(i)^2); 
    %   Laplace Approx 
    for t = 1:T
        ht = betat; 
        gra = 1;
        while abs(gra)>1e-06
            gra = 0.5 - .5*exp(-ht)*s2(t,i) + (1/gamma2)*(ht-betat);
            hes = 1/gamma2 + .5*exp(-ht)*s2(t,i);
            ht = ht - gra/hes;
        end
        hfil(t,i) = ht; 
        % betats(t,i) = betat; % gamma2s(t,i) = gamma2; % sigma2s(t,i) = 1./hes;
        llh_t(t,i) = 0.5*log(2*pi) - (.5*log(2*pi*exp(ht))+.5*exp(-ht)*s2(t,i)+.5*log(2*pi*gamma2)+ .5*(1/gamma2)*(ht-betat)^2) - 0.5.*log(abs(hes));
        % Update predictive
        betat =phi(i)*ht;
        gamma2 = phi(i)^2./hes + s(i);
    end  
end 
llh = -log(abs(det(B))) + sum(llh_t,2);     
    
end

 
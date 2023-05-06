function llh = Fh1c( param_vec,  y, p, S_b,sumdum)
% Log Likelihood Based on Laplace Approximation at h=h_lap
%% Calculate negative LLH
Y = y(p+1:end,:);
[T,K] = size(Y);
yl = lagmatrix(y,1:p);
X = [ones(T,1),yl(p+1:end,:)];
nA = K*(K*p+1);
vAlp = param_vec(1:nA);
vB = param_vec(nA+1:nA+K^2);
phi = param_vec(nA+K^2+1:nA+K^2+K);
s = param_vec(nA+K^2+K+1:end);
A = reshape(vAlp,K,(K*p+1))'; % A prime!
u = Y-X*A; 
B = reshape(S_b'*vB,K,K);
e = (B\u')'; 
s2 = e.^2+ std(e)./20; 
muh = -.5* s./(1-phi.^2);  

%% Filter algorithm to evaluate log likelihood for each shock
llh_t = zeros(T,K);
hfil = zeros(T,K); 
%gamma2s = zeros(T,K); betats = zeros(T,K); sigma2s = zeros(T,K);
for i = 1:K 
    betat = muh(i);
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
        betat = muh(i) + phi(i)*(ht-muh(i));
        gamma2 = phi(i)^2./hes + s(i);
    end  
end 
llh = -log(abs(det(B))) + sum(llh_t,2);  
if sumdum == 1
    llh = sum(llh);
end
 
end
















% % 
% % 
% % 
% % %% Newton Raphson  
% % llh_t = zeros(T,K);
% % hfil = zeros(T,K); 
% % for i = 1:K 
% %  % For t=1 we have: 
% % h1 = h_lap(1,i);
% % gra = 1; 
% % while abs(gra)>1e-06
% %     gra = 0.5 - .5*exp(-h1)*s2(1,i)  + (1-phi(i)^2)/s(i)*(h1-muh(i));
% %     hes = (1-phi(i)^2)/s(i) + .5*exp(-h1)*s2(1,i) ;
% %     h1 = h1 - gra/hes;  
% % end
% % hfil(1,i) = h1; 
% % llh_t(1,i) = 0.5*log(2*pi) - (.5*log(2*pi*exp(h1))+.5*exp(-h1)*s2(1,i) +.5*log(2*pi*s(i)/(1-phi(i)^2))+ .5*(1-phi(i)^2)/s(i)*(h1-muh(i))^2) - 0.5.*log(abs(hes));
% % sig2m1 = 1./hes; 
% % htm1 = h1;
% % for t = 2:T 
% %    betat = muh(i) + phi(i)*(htm1-muh(i));
% %    gamma2 = phi(i)^2*sig2m1 + s(i); 
% %    gra = 1;
% %    ht = betat;
% %    while abs(gra)>1e-06
% %        gra = 0.5 - .5*exp(-ht)*s2(t,i) + (1/gamma2)*(ht-betat);
% %        hes = 1/gamma2 + .5*exp(-ht)*s2(t,i);
% %        ht = ht - gra/hes;
% %    end
% %    hfil(t,i) = ht;
% %    llh_t(t,i) = 0.5*log(2*pi) - (.5*log(2*pi*exp(ht))+.5*exp(-ht)*s2(t,i)+.5*log(2*pi*gamma2)+ .5*(1/gamma2)*(ht-betat)^2) - 0.5.*log(abs(hes));
% %    sig2m1 = 1./hes;  
% %    htm1 = ht;
% % end 
% % end
% % llh = -log(abs(det(B))) + sum(llh_t,2); 
% % 
% % if sumdum == 1
% %     llh = sum(llh);
% % end

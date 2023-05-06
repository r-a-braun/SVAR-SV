function [llh, intlikestd, pval] = llh_IS( u, B, M, phi, s, H)
% computes likelihood based on importance sampling
epsilon = (B\u')';   
[T,K] = size(epsilon); 
idx_r = find(isnan(phi)==0);
idx_nr = find(isnan(phi));
r = length(idx_r);
[llh_heterosk, intlikestd, pval] = intlike_svsvar(epsilon(:,idx_r),s(idx_r),phi(idx_r),H(:,idx_r),M);  
llh_homosk = - T*(K-r)/2*log(2*pi) - 1/2* (sum(sum(epsilon(:,idx_nr).^2)));
llh = - T*log(abs(det(B))) + llh_heterosk + llh_homosk ; 
end

% This function relies on code of Chan and Eisenstat (2017)
% See:
% Chan, J.C.C. and Eisenstat, E. (2017). Bayesian model comparison for
% time-varying parameter VARs with stochastic volatility, Journal of
% Applied Econometrics, forthcoming.

function [intlike, intlikestd, pval_finitevariance] = intlike_svsvar(struct_err,s,phi,hinit,R)
[T,r]=size(struct_err);
e = vec(struct_err'); 
% obtain the proposal density
Hh = speye(T*r) - sparse(r+1:T*r,1:(T-1)*r,repmat(phi,T-1,1),T*r,T*r);
HinvSH_h = Hh'*sparse(1:T*r,1:T*r, [(1-phi.^2)./s;repmat(1./s,T-1,1)])*Hh; 
s2 = vec((struct_err.^2 + std(struct_err)./1000)'); 
ht = vec(hinit'); 
%% Newton Raphson
A = kron(ones(1,T)/T,eye(r)); % linear constraint: Ax = 0 
e_h = inf;
while e_h >  1e-06
    einvhts2 = exp(-ht).*s2;
    gh = -HinvSH_h*(ht) -.5*(1-einvhts2);
    Gh = -HinvSH_h -.5*sparse(1:T*r,1:T*r,einvhts2);             
    newht = ht - Gh\gh; 
    W = -Gh\A';
    newht = newht - W*((A*W)\(A*newht)); 
    e_h = max(abs(newht-ht));
    ht = newht;   
end      
Kh = -Gh; 
CKh = chol(Kh,'lower');

%% Draw from the IS
bigh = repmat(ht,1,R) + CKh'\randn(T*r,R);
WAW = W/(A*W);
bigh = bigh - WAW*(A*bigh); 

%% IS estimator of the likelihood:  
cons_IS = -T*r/2*log(2*pi) + sum(log(diag(CKh))) + (- 0.5*log(det(A*A'))) - ( -r/2*log(2*pi) - sum(log(diag(chol(A*W,'lower')))));  
cons_LLH = -T*r/2*log(2*pi);
cons_prior = -T*r/2*log(2*pi) - 0.5*log(det(A*A')) + r/2*log(2*pi);  
logIS = cons_IS   -.5*sum(((bigh-ht)'*Kh)'.*(bigh-ht)) ;
% Likelihood  
logLike =  cons_LLH - .5*sum(bigh) - .5*sum(e.^2.*exp(-bigh)); 
% Prior
logPrior = cons_prior - sum(.5*(T*log(s) - log(1-phi.^2)))-.5*sum(((bigh)'*HinvSH_h)'.*(bigh))...
    + sum(log(diag(chol((A*(HinvSH_h\A')),'lower')))); 
store_llike = (logLike + logPrior - logIS)'; 
maxllike = max(store_llike);
intlike = log(mean(exp(store_llike-maxllike))) + maxllike;
[~,pval_finitevariance] = FiniteCheck(exp(store_llike-maxllike),.5,.05); 
R = floor(R/20)*20; % in case that its not a multiple of 20
store_llike = store_llike(1:R);
shortw = reshape(store_llike,R/20,20);
maxw = max(shortw);
bigml = log(mean(exp(shortw-repmat(maxw,R/20,1)),1)) + maxw;
intlikestd = std(bigml)/sqrt(20); 
end
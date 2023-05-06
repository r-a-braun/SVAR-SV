function [ Param_est, loglike, count, VarTheta, SAnalysis] = SVAR_SV_em_approx( y, specs)
% Bertsche, Braun (2020): Identification of Structural Vector
% Autoregressions by Stochastic Volatility
% Main File that implements 
% 1) Runs EM-1 Algorithm  
% 3) Computes Loglikelihood values
% 4) Computes estimate of asymptotic covariance matrix based on Oakes method
% 5) Computes bootstrap-type confidence intervals of IRFs and FEVDs based on 2999 simulations from the asymptoric
% distribution

%% Model:
% SV-SVAR model restricte s.t. A*h = mu_i = 0
% h_ti = mu_i+ phi_i* h_{t-1,i} + w_{ti}      w_{ti} ~ N(0,s_i);
% y_t = A x_t + B exp(h_{t})*eta_{t}           eta_{t}  ~ N(0,I_n);

%% Read Specifications:
p = specs.p;
Param_init = specs.param_init; 
M = specs.M; 
max_it = specs.max_it;
r = specs.r; 
tol = specs.tol;
Y = y(p+1:end,:);
[T,K] = size(Y); 
S_b = eye(K^2); 
Restr = ones(K,K); Restr(r+1:end,r+1:end) =tril(ones(K-r,K-r));
S_b(:,vec(Restr)==0)=[]; 
m = K*(K*p+1);
scale = specs.scale;
yl = lagmatrix(y,1:p);
X = [ones(T,1),yl(p+1:end,:)]; 
if isempty(Param_init) % Generate some Starting Values: 
    [Ahat_ls, Sigmahat, u_it] = VARmlike(y,p); 
    [Qran,~] = qr(randn(K,K)); 
    B_it = chol(Sigmahat)'*Qran; 
    B_it = B_it*RotateIdent(B_it,r); 
    vB_it = pinv(S_b)*B_it(:);
    B_it = reshape(S_b*vB_it,K,K);
    A_it = Ahat_ls';
    e_it = (B_it\u_it')';
    phi_it = NaN(K,1); s_it=NaN(K,1);
    phi_it(1:r) = repmat(0.95,r,1);
    s_it(1:r) = repmat(0.02,r,1);
    Ht = zeros(T,K);
else
    A_it = Param_init.A;
    vA_it = A_it(:);
    B_it = Param_init.B;
    vB_it = pinv(S_b)*B_it(:);
    B_it = reshape(S_b*vB_it,K,K);
    phi_it  = NaN(K,1); s_it=NaN(K,1);
    phi_it(1:r) = Param_init.phi(1:r) ;
    s_it(1:r) = Param_init.s(1:r);
    A_it = reshape(vA_it,m/K,K);
    u_it = Y - X*A_it;
    e_it = (B_it\u_it')';
    Ht = Param_init.H;
end

%% Preliminaries
options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off','MaxFunEvals', 10000 );   
Param_est.X = X;
Param_est.Y = Y;
Z = kron(speye(K),X);
ytil = Y(:);  

%% EM Algorithm
llh_diff = inf;
count = 0;  
VarHt = zeros(T,K); CovHtHtm1 = zeros(T-1,K); iV = ones(T,K); 
idx_r = find(isnan(phi_it)==0)';
while and(llh_diff>tol, count<max_it)
    count = count +1 ;  
    %% E Step and M step for phi and s
    dQ1s = zeros(K,1);
    for j=idx_r
        [ Ht(:,j), VarHt(:,j), CovHtHtm1(:,j), iV(:,j) ] = e_step_approx( e_it(:,j), phi_it(j), s_it(j), Ht(:,j),scale);
        [ phi_it(j), s_it(j), dQ1s(j) ] = m_step_phis(phi_it(j), s_it(j), Ht(:,j), VarHt(:,j), CovHtHtm1(:,j));
    end 
    dQ1 = sum(dQ1s);  
    %% M step for A and B
    [ B_it, A_it, u_it, e_it, dQ2 ] = m_step_AB(B_it, A_it, iV, Y, X,ytil, Z, S_b,options);
    llh_diff = dQ1 + dQ2;  
   % disp(llh_diff) 
end    
%% Order the shocks as to maximize the trace
[ pvec, signs ] = perm_vec( B_it(1:r,1:r) );
B_it(:,1:r) = B_it(:,pvec)./signs;
Ht(:,1:r) = Ht(:,pvec); 
phi_it(1:r) = phi_it(pvec);
s_it(1:r) = s_it(pvec);
e_it(:,1:r) = e_it(:,pvec)./signs;
VarHt(:,1:r) = VarHt(:,pvec); 
iV(:,1:r) = iV(:,pvec); 
CovHtHtm1(:,1:r) = CovHtHtm1(:,pvec);
idx_r = find(isnan(phi_it)==0)';


%% Save output
Param_est.A = A_it; 
Param_est.u = u_it;
Param_est.phi = phi_it; 
Param_est.B = B_it;
Param_est.s = s_it; 
Param_est.H = Ht;
Param_est.e = e_it;
Param_est.VarH = VarHt;
  
%% EVALUATE likelihood if requested
if nargout > 1  
    [llh, intlikestd,pval_fv] = llh_IS( u_it, B_it, M, phi_it, s_it, Ht);
    loglike.llh = llh;
    loglike.pval_fv = pval_fv;
    loglike.llh_sdv = intlikestd;
    n_param = length(s_it) + length(phi_it) + length(vB_it) + length(vec(A_it)); 
    loglike.aic = -2*llh + 2*n_param; 
    loglike.bic = -2*llh + log(T)*n_param;
end  

%% IF REQUESTED: COMPUTE STANDARD DEVIATIONS by numerical differentiation
if nargout>3  
    nA = K*(K*p+1); 
    paramhat = [phi_it(idx_r);s_it(idx_r);vec(A_it');S_b'*vec(B_it)];
    Qd2 = Expected_Hess(Param_est,Y,X,Ht,VarHt,CovHtHtm1,iV,S_b);
    objfun_g = @(x)Expected_Grad(x,Param_est,Y,X,p,Ht,S_b,scale);  
    Qd1d2 = Gradp(objfun_g,paramhat);
    HessianOakesFormula = Qd2 + Qd1d2;%(Qd1d2+Qd1d2')/2;
    Var_hat = inv(-HessianOakesFormula);   
    [~,pchol]=chol(Var_hat,'lower');
    if pchol ~= 0 % Getting Qd1d2 pos. sem. can be a numerical pain if many parameters (try increasing M first)
        VarTheta.cov =  (Var_hat + Var_hat')/2; 
    else
        VarTheta.cov = Var_hat; 
    end 
    VarTheta.phi =  Var_hat(1:r, 1:r);
    VarTheta.s   =  Var_hat(r+1:2*r, r+1:2*r);
    VarTheta.vA  =  Var_hat(2*r+1:2*r+nA, 2*r+1:2*r+nA);
    VarTheta.vB  =  Var_hat(2*r+nA+1:end ,2*r+nA+1:end);  
    VarTheta.vAvB = Var_hat(2*r+1:2*r+nA, 2*r+nA+1:end);
end
if nargout > 4
    hor = specs.hor;   
    %% Compute IRFs and FEVD based on parametric "bootstrap"  (1999 simulations) 
    V = exp( 0.5*s_it./(1-phi_it.^2) );
    Bimp = B_it*diag( sqrt(V) );
    [AA, J] = companion(A_it', 1);  
    SVMA = zeros(K,K,hor+1);
    THETA = zeros(K,K,hor+1);
    FEVD = zeros(K,K,hor+1);
    RS = eye(K);
    if specs.unitshocks == 1
        RS = diag(1./diag(Bimp)); 
    end 
    for ii = 0:hor
        SVMA(:,:,ii+1) = J*(AA^ii)*J'*Bimp;
        FEVD(:,:,ii+1) = sum(SVMA(:,:,1:ii+1).^2,3)./sum(sum(SVMA(:,:,1:ii+1).^2,3),2);
        THETA(:,:,ii+1) = SVMA(:,:,ii+1)*RS;
    end 
    PSI = cumsum(THETA,3);  
    SAnalysis.irfs = reshape(THETA,K^2,hor+1)';
    SAnalysis.irfs_cs = reshape(PSI,K^2,hor+1)';
    SAnalysis.fevd = reshape(FEVD,K^2,hor+1)';
    
    % Start bootstrap iterations
    nboot = 2999; 
    irf_save  = zeros(hor+1,K^2,nboot);
    cirf_save = zeros(hor+1,K^2,nboot);
    FEVD_save = zeros(hor+1,K^2,nboot);
    for i = 1:nboot
        param_i = paramhat + chol(VarTheta.cov,'lower')*randn(length(paramhat),1); % Draw a model from the asymptotic distribution
        phii = param_i(1:r);
        si = param_i(r+1:2*r);
        vAlpi = param_i(2*r+1:2*r+nA);
        vBi = param_i(2*r+nA+1:end);
        Ai = reshape(vAlpi,K,(K*p+1))'; 
        var_shocks = [exp(0.5*si./(1-phii.^2));ones(K-r,1)];
        V12 = diag(sqrt(var_shocks));
        Bimp = reshape(S_b*vBi,K,K)*V12;
        [AA, J] = companion(Ai', 1);
        SVMA = zeros(K,K,hor+1); THETA = zeros(K,K,hor+1); FEVD = zeros(K,K,hor+1);
        RS = eye(K);
        if specs.unitshocks == 1
            RS = diag(1./diag(Bimp));
        end
        for ii = 0:hor
            SVMA(:,:,ii+1) = J*(AA^ii)*J'*Bimp;
            FEVD(:,:,ii+1) = sum(SVMA(:,:,1:ii+1).^2,3)./sum(sum(SVMA(:,:,1:ii+1).^2,3),2);
            THETA(:,:,ii+1) = SVMA(:,:,ii+1)*RS;
        end
        PSI = cumsum(THETA,3);
        irf_save(:,:,i) = reshape(THETA,K^2,hor+1)';
        cirf_save(:,:,i) = reshape(PSI,K^2,hor+1)';
        FEVD_save(:,:,i) = reshape(FEVD,K^2,hor+1)';
    end
            
 
    for ii=1:size(specs.alpha,2) 
        SAnalysis.IRF_qu(:,:,ii) = SAnalysis.irfs + quantile(irf_save-SAnalysis.irfs ,1-specs.alpha(ii)/2,3);
        SAnalysis.IRF_ql(:,:,ii) = SAnalysis.irfs + quantile(irf_save-SAnalysis.irfs , specs.alpha(ii)/2,3);
        SAnalysis.IRF_cs_qu(:,:,ii) = SAnalysis.irfs_cs + quantile(cirf_save-SAnalysis.irfs_cs ,1-specs.alpha(ii)/2,3);
        SAnalysis.IRF_cs_ql(:,:,ii) = SAnalysis.irfs_cs + quantile(cirf_save-SAnalysis.irfs_cs , specs.alpha(ii)/2,3);
        SAnalysis.fevd_qu(:,:,ii) = SAnalysis.fevd + quantile(FEVD_save-SAnalysis.fevd ,1-specs.alpha(ii)/2,3);
        SAnalysis.fevd_ql(:,:,ii) = SAnalysis.fevd + quantile(FEVD_save-SAnalysis.fevd , specs.alpha(ii)/2,3);
        
    end  
end

end







    
    
    
%     
%     
%     dQ1 = 0;
%     for j=1:K
%         h_ks = H_i(:,j); 
%         P_ks = Var_H(:,j);  
%         Ptm1t_s = Cov_Ht1(:,j);
%         Sxx = sum(P_ks(1:end-1) + h_ks(1:end-1).^2);
%         Sxy = sum( Ptm1t_s+h_ks(1:end-1).*h_ks(2:end) );
%         Syy = sum(P_ks(2:end)+h_ks(2:end).^2);
%         beta_old = phi_it(j)';
%         I2old = -0.5*( (T-1)*log(s_it(j)) +  1/s_it(j)*(Syy - 2*Sxy'*beta_old  + beta_old'*Sxx*beta_old));
%         beta_max = Sxx\Sxy;
%         s_it(j) = 1/(T-1)*(Syy - 2*Sxy'*beta_max  + beta_max'*Sxx*beta_max); 
%         phi_it(j) = beta_max; 
%         I2new = (-0.5*((T-1)*log(s_it(j)) + (T-1)));
%         dQ1 = dQ1 + (I2new-I2old);
%     end  


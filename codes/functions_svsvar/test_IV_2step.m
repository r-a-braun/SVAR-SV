function [ output ] = test_IV_2step( y,   input )
%% Input arguments:

param_est = input.param_est;
e = param_est.e;
B = param_est.B; 
A = param_est.A; 
s = param_est.s;
phi = param_est.phi; 
p = input.specs.p;
m = input.m; 
nm = size(m,2);
theta1step = [vec(A'); vec(B)];
Y = y(p+1:end,:);
[T,K] = size(Y);   
yl = lagmatrix(y,1:p);
X = [ones(T,1),yl(p+1:end,:)];  
df1dth1 = Gradp(@(x)prederr_llh(Y, X ,p, phi, s, x ),  theta1step )' ; 
V1 = T*input.VarTheta.cov(2*K+1:end,2*K+1:end); 
for j = 1:nm
    idxj = isnan(m(:,j))==0;
    ej = e(idxj,:);
    mj = m(idxj,j);
    Tj = length(ej); 
    phih = inv(ej'*ej)*ej'*mj; 
    sige2h = (mj - ej*phih)'*(mj - ej*phih)/(Tj);
    V2 = sige2h*inv(ej'*ej/(Tj));
    df2dth1 =  Gradp(@(x)llh2ndstep( Y(idxj,:) , X(idxj,:), p, mj , sige2h , x, phih ),theta1step)';
    df2dth2 = Gradp(@(x)llh2ndstep( Y(idxj,:) , X(idxj,:), p, mj, sige2h , theta1step, x ),phih)'; 
    C = 1/(Tj)*df2dth2*df2dth1';
    R = 1/(Tj)*df2dth2*df1dth1(:,idxj)';
    V2_2step = (V2 + V2*( C*V1*C' - R*V1*C' - C*V1*R')*V2);  
    [maxphi,idxmax]=max(abs(phih));
    rest = setdiff(1:K,idxmax);  
    output.phi(:,j) = phih;
    output.se(:,j) = sqrt(diag(V2_2step)/T);
    output.corr(:,j) = corr(ej,mj);
    output.maxcorr(j) = corr(ej(:,idxmax),mj);
    output.nonrelevance(j) = (1-tcdf(maxphi./sqrt(V2_2step(idxmax,idxmax)/T),Tj))*2;
    output.exogeneity(j) = 1-chi2cdf(phih(rest)'*inv(V2_2step(rest,rest)/T)*phih(rest),length(rest));
    output.indices(j) = idxmax;
end
  
end

function llh = llh2ndstep( Y, X, p, m, sig2 , theta1, theta2 )
%LLH1 Summary of this function goes here
%   Detailed explanation goes here
[T,K]=size(Y); 
nA = K*(K*p+1);   
vAlp = theta1(1:nA);
vB = theta1(nA+1:end);  
A = reshape(vAlp,K,(K*p+1))'; 
B = reshape(vB,K,K);   
E = (B\(Y - X*A)')';
phi = theta2(1:K);  
e2 = m - E*phi; 
llh = -1/2*log(det(sig2)) - 0.5*(e2.^2/sig2);
end


 


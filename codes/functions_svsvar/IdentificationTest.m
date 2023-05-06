function TestResults = IdentificationTest(y,ic,p,hor)
% Bertsche, Braun (2020): Identification of Structural Vector
% Autoregressions by Stochastic Volatility
%
% function to evaluate test statistics and give p-values (Lanne & Saikkonen
% (2007))
%
%% Input:
% - y (T+p x K): raw time series vector
% - p (1 x 1): number of autoregressive lags
% - ic (1 x 1): 1 if intercept, 0 if not
% - A (m x K): intercept and AR matrices
% - B (K x K): impact matrix
% - hor (1 x 1): considered horizon for autocorrelations
%
%% Output:
% - TestResults (structure)

input.y = y; input.p = p; input.ic = ic;
[Y,X] = VAR_Trans(y,p,ic); %data transformation
[T,K] = size(Y); %time series length and number of variables
M = VECH(K);
Q1 = zeros(K,1); Q2 = Q1;
dofQ1 = hor.*ones(K,1); dofQ2 = zeros(K,1);
for r = 0:K-1
    M = VECH(K-r);
    input.r = r;
    Parameters = Est_SV_EM1(input);
    u = Y-X*Parameters.A; %reduced form residuals
    eps = (Parameters.B\u')'; %structural shocks
    eps2 = eps(:,r+1:end); %consider homoskedastic shocks
    xi_t = sum(eps2.^2,2);
    xi = xi_t-mean(xi_t);
    nu_t = zeros((K-r)*(K-r+1)/2,T);
    for t = 1:T
        temp = eps2(t,:)'*eps2(t,:);
        nu_t(:,t) = M*temp(:);
    end
    nu = nu_t-mean(nu_t,2); 
    gamma0 = mean(xi.^2);
    Gamma0 = (nu*nu')./T;
    iG0 = inv(Gamma0);
    gammah = zeros(hor,1);
    for i = 1:hor
        gammah(i) = sum(xi(i+1:end).*xi(1:end-i))./T;
        Gammah = (nu(:,i+1:end)*nu(:,1:end-i)')./T;
        Q2(r+1) = Q2(r+1)+trace(Gammah'*iG0*Gammah*iG0);
    end
    Q2(r+1) = T*Q2(r+1);
    Q1(r+1) = T.*sum((gammah./gamma0).^2);
    dofQ2(r+1) = hor/4.*(K-r)^2*(K-r+1)^2;
end
TestResults.Q1.stat = Q1;
TestResults.Q1.dof = dofQ1;
TestResults.Q1.pv = 1-chi2cdf(Q1,dofQ1);

TestResults.Q2.stat = Q2;
TestResults.Q2.dof = dofQ2;
TestResults.Q2.pv = 1-chi2cdf(Q2,dofQ2);
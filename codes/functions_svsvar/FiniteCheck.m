function [ind,pv] = FiniteCheck(weights,tv,alpha)
% Bertsche, Braun (2020): Identification of Structural Vector
% Autoregressions by Stochastic Volatility
%
% function to test finiteness of importance weights (Koopman, Shephard & Creal (2009))
%
%% Input:
% - weights (R x r): importance weights
% - tv (1 x 1): percentile for threshold of check for finiteness of
% importance weights
% - alpha (1 x 1): confidence level
%
%% Output:
% - ind (1 x 1): 1 if weights are finite, 0 else

%% Specify optimization options
options =optimoptions('fmincon','Algorithm','interior-point','Display','off',...
        'MaxFunctionEvaluations',1e+05,'StepTolerance',1e-15);

%% Side conditions of optimization
A = -eye(2); b = [-1e-03;-0.5];
[~,r] = size(weights); %number of heteroskedastic shocks
pv = zeros(r,1); %empty matrix for p-values
for i = 1:r
    thresh = quantile(weights(:,i),tv); %threshold value: quantile of weights
    id = weights(:,i)>thresh; %consider weights above threshold
    z = weights(id,i)-thresh; %consider weights above threshold
    %% Unconstrained ML estimation
    [~,L1] = fmincon(@(theta)neg_log_GPD(theta,z),[1;0.8],A,b,[],[],[],[],[],options);
    %% ML estimation conditional on xi=0.5
    [~,L0] = fmincon(@(theta)neg_log_GPD(theta,z),[1;0.5],A,b,[0,1],0.5,[],[],[],options);
    %% LR test statistic and p-value
    LR = 2.*(L1-L0);
    pv(i) = .5-0.5.*chi2cdf(LR,1);
end  
ind = prod(pv>alpha); %if all p-values are above alpha --> set ind =1, 0 else



end



function val = neg_log_GPD(theta,z)
% Bertsche, Braun (2018): Identification of Structural Vector
% Autoregressions by Stochastic Volatility
%
% function to evaluate log pdf of generalized pareto distribution
%
%% Input:
% - theta (2 x 1): vector of distribution parameters
% - z (tv.*R x 1): observations
%
%% Output:
% - val (1 x 1): log pdf value

N = length(z); %number of observations
beta = theta(1); xi = theta(2); %distribution's parameters
val = N*log(beta)+(1+1./xi).*sum(log(1+xi.*z./beta)); %log pdf value
end

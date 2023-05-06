function y_sim = Simulate_SV(A,B,phi,s,T,bi,rep,mink)
% Bertsche, Braun (2020): Identification of Structural Vector
% Autoregressions by Stochastic Volatility
%
% function to simulate time series with SV structural errors
%
%% Input:
% - A (K x Kp+ic): intercept and VAR matrix
% - B (K x K): impact matrix
% - phi (1 x 1): AR parameter SV equation
% - s (1 x 1): error term variance SV equation
% - T (1 x 1): time series length
% - bi (1 x 1): burn-in
% - rep (1 x 1): no. of replications
% - mink (1 x 1): minimum kurtosis
%
%% Output:
% - y_sim (K x T x rep): simulated time series

%% Figure out size and Parameters
[K,Kpic] = size(A);
p = floor(Kpic./K); ic = Kpic-K*p;
if ic == 1
    IC = A(:,1);
else
    IC = zeros(K,1);
end
Avar = A(:,ic+1:end);

%% Simulate error terms
e_M = randn(K,T+bi,rep);
muh = 0; 
Hphi = speye(T+bi) - sparse(2:T+bi,1:(T+bi-1),phi*ones(1,T+bi-1),T+bi,T+bi);
HiSH = Hphi'*spdiags([(1-phi^2)/s; 1/s*ones(T+bi-1,1)],0,T+bi,T+bi)*Hphi;
deltah = Hphi\([muh;muh*(1-phi)*ones(T+bi-1,1)]); 
A_r = ones(1,T+bi)/(T+bi);
W = HiSH\A_r';   
WAW = W/(A_r*W);  
cholKh = chol(HiSH,'lower'); 
H = zeros(T+bi,rep,K); idx = (1:rep); bigh = zeros(T+bi,rep);
while isempty(idx) == 0
    for k = 1:K
        bigh(:,idx) = repmat(deltah,1,length(idx)) + cholKh'\randn(T+bi,length(idx));
        H(:,idx,k) = bigh(:,idx) - WAW*(A_r*bigh(:,idx) - muh);   
    end
    for t= 2:T+bi
        H_per = reshape(permute(H(t,idx,:),[3,2,1]),K,1,length(idx));
        e_M(:,t,idx) = exp(0.5*H_per).*randn([K,1,length(idx)]); 
    end
    temp = zeros(1,rep);
    for k = 1:K
        temp1 = squeeze(e_M(k,bi+1:end,:));
        temp2 = kurtosis(temp1);
        temp = temp+(temp2<mink);
    end
    idx = find(temp~=0);
end

eps_m = reshape(e_M,K,rep*(T+bi));  
u_M = (B*eps_m)'; 
u_M = reshape(u_M',K,T+bi,rep);  
if rep > 1
    u_M = u_M- mean(u_M,3);
end

%% Simulate time series
y_sim = zeros(K,T+bi,rep); 
y_sim(:,1:p,:)= repmat(IC*ones(1,p),1,1,rep) + u_M(:,1:p,:); 
v_3d = repmat(IC,1,1,rep); 
for t = p+1:(T+bi)
    y_sim(:,t,:)= v_3d + reshape(Avar* reshape(y_sim(:,t-1:-1:t-p,:),K*p,rep),K,1,rep)  + u_M(:,t,:);
end  
y_sim = y_sim(:,bi+1:end,:);
end
function llh = F2( y, X, p, m, sig2 , theta1, theta2 )
%LLH1 Summary of this function goes here
%   Detailed explanation goes here
[T,K]=size(y);
Alpha = reshape(theta1(1:K*(K*p+1)),K,K*p+1);
B = reshape(theta1(K*(K*p+1)+1:end),K,K); 
E = (B\(y - X*Alpha')')';
phi = theta2(1:K);  
e2 = m - E*phi; 
llh = -1/2*log(det(sig2)) - 0.5*(e2.^2/sig2);
end


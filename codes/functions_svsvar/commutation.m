function K = commutation(m, n)
% Constructs Commutation Matrix K_mn as defined in Lütkepohl p. 663
a = zeros(m*n, m);  
for i = 1:m 
    j = (i - 1) * n + 1; 
    a(j, i) = 1; 
end 
K = a; 
for i = 2:n 
    K = [K, [zeros(i-1, m); a(1:end-i+1, :)]]; 
end
K = sparse(K);

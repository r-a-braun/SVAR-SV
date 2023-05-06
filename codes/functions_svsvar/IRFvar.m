function IRFs = IRFvar( A, B,hor,inc)
% Computes Impulse response functions up to horizon hor
K = size(B,2); 
[ AA, J , ~] = companion( A , inc);
[ pvec, signs ] = perm_vec( B);
B = B(:,pvec)./signs;
Bunit = B*diag(1./diag(B));
Phis = zeros(K,K,hor+1);
Thetas = zeros(K,K,hor+1); 
Phis(:,:,1) = eye(K);
Thetas(:,:,1) = Bunit; 
Apoweri = AA;
for i = 1:hor
    Phis(:,:,i+1) = J*Apoweri*J';
    Thetas(:,:,i+1) = Phis(:,:,i+1)*Bunit;
    Apoweri = Apoweri * AA;
end 
IRFs.irfs = reshape(Thetas,K^2,hor+1)';
IRFs.irfs_cs = cumsum(IRFs.irfs);

    
end


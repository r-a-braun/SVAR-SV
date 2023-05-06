function [ B, A, U, E, dQ2 ] = m_step_AB(B, A, iVar_h, Y, X, ytil, Z, S_b, options) 
%% 2) Update A (closed form) & B (numerical)
[T,K]=size(Y);
vB = S_b'*vec(B);
invSig = sparse(1:T*K,1:T*K, iVar_h(:) ) ;
nQ2_old = min_logpost_B( vB, iVar_h, Y-X*A, S_b);
% Update A
iB = inv(B);
iBs = kron(iB,speye(T));
iSkI = iBs'*invSig*iBs;
iV_p = Z'*iSkI*Z ;
ahat = iV_p\(Z'*iSkI*ytil );
A = reshape(ahat,size(A,1),size(A,2));
U = Y-X*A;
% Update B 
[vB,nQ2_new] = fminunc(@(vB)min_logpost_B(vB, iVar_h, U , S_b),vB,options);  
B = reshape(S_b*vB,K,K);
E = (B\U')';
dQ2 = - nQ2_new  - (-nQ2_old) ;
end
function [ nllh ]  = min_logpost_B( vB,  iVar_h,  U,  S_b) 
[T,K] = size(U); 
B = reshape(S_b*vB,K,K);  
E = (B\U')'; 
nllh = T*log(abs(det(B)))+0.5*sum(sum(E.^2.*iVar_h));  
end
function [h_ks, P_ks, Ptm1t_s, iVar_h] = e_step_approx(struct_err, phi, s, hinit ,scale) 
% E STEP APPROXIMATION
[T,r]=size(struct_err);
e = vec(struct_err'); 
Hh = speye(T*r) - sparse(r+1:T*r,1:(T-1)*r,repmat(phi,T-1,1),T*r,T*r); 
HinvSH_h = Hh'*sparse(1:T*r,1:T*r, [(1-phi.^2)./s;repmat(1./s,T-1,1)])*Hh;
s2 = vec((struct_err.^2 + std(struct_err)./scale)'); 
e_h = 1; 
ht = vec(hinit'); 
count = 0;
%% Newton Raphson
A = kron(ones(1,T)/T,eye(r)); % linear constraint: Ax = 0 
while e_h >  1e-07
    einvhts2 = exp(-ht).*s2;
    gh = -HinvSH_h*(ht) -.5*(1-einvhts2);
    Gh = -HinvSH_h -.5*sparse(1:T*r,1:T*r,einvhts2);             
    newht = ht - Gh\gh; 
    W = -Gh\A';
    newht = newht - W*((A*W)\(A*newht)); 
    e_h = max(abs(newht-ht));
    ht = newht;  
    count = count + 1;
end      
Kh = -Gh;  
[Ch, ~, Q] = chol(Kh,'lower');  
P = Q';
d = diag(Ch) ; 
Ch = tril(Ch / diag (d), -1) ;
d = d.^2 ;
d = full(d) ;
[Z , ~] = sparseinv_mex(Ch, d, Ch, spones(P*Kh*Q)); %Takayashis equations  
Sigma = Q*Z*P ;   % Inverse only computed for relevant elements ( variance and lag 1 covariance)
% spy(Sigma)
%% Correct for linear constraints
C = W*((A*W)\W'); 
%% Compute Moments
h_ks = reshape(ht,r,T)'; % - mean(ht);
P_ks = reshape(diag(Sigma) - diag(C),r,T)';
Ptm1t_s = reshape( diag(Sigma,-r) - diag(C,-r) ,r,T-1)';    
iVar_h = exp(- h_ks + 0.5*P_ks );  


end
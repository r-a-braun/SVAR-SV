function [ H_i,Var_H,Cov_Htm1,iVar_h ] = e_step_IS( e, phi, s, hinit, CRN ,scale )
%E_STEP_IS Summary of this function goes here
%   Detailed explanation goes here

[T,nrep]=size(CRN);  
r = size(e,2);
Hh = speye(T*r) - sparse(r+1:T*r,1:(T-1)*r,repmat(phi,T-1,1),T*r,T*r); 
HinvSH_h = Hh'*sparse(1:T*r,1:T*r, [(1-phi.^2)./s;repmat(1./s,T-1,1)])*Hh;
s2 = (e).^2  + std(e)./scale;
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
CKh = chol(Kh,'lower'); 
WAW = W/(A*W); 
vHsim = ht + CKh'\CRN;
vHsim = vHsim - WAW*(A*vHsim); 
%% Importance Density estimator
store_llike = - .5*sum(vHsim) - .5*sum(s2.*exp(-vHsim))  -.5*sum(((vHsim)'*HinvSH_h)'.*(vHsim))...
    -(-.5*sum(((vHsim-ht)'*Kh)'.*(vHsim-ht))); 
maxllike = max(store_llike);
weights = exp(store_llike-maxllike);
weights = weights./sum(weights);  
H_i = vHsim*weights';
Var_H = (vHsim-H_i).^2*weights';
Cov_Htm1 = ((vHsim(1:end-1,:)-H_i(1:end-1,:)).*(vHsim(2:end,:)-H_i(2:end,:)))*weights';
iVar_h = exp(-vHsim)*weights';
end


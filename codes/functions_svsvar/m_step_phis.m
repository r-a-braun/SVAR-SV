function [ phi, s, dQ1 ] = m_step_phis(phi, s, H_i, Var_H, Cov_Ht1 )
%M_STEP_PHIS Summary of this function goes here
%   Detailed explanation goes here
[T,r]=size(H_i);
vSxx = sum(Var_H(1:end-1,:) + H_i(1:end-1,:).^2)';
vSxy = sum(Cov_Ht1 + H_i(1:end-1,:).*H_i(2:end,:) )';
vSyy = sum(Var_H(2:end,:)+H_i(2:end,:).^2)';
Eh1sq = (Var_H(1,:) + H_i(1,:).^2)';
[I2old,gg,hh] = Grad_Hess_eprior([phi;s],vSyy,vSxx,vSxy,Eh1sq,T,r);
phi = vSxy./vSxx; 
if sum(phi>.99)>0
    phi=.99;
end
s = 1/(T-1).*(vSyy - 2.*vSxy.*phi  + phi.^2.*vSxx);
ps = [phi;s];
dd = inf;
count1 = 0; stepsize = 1;
while dd > 1e-05 && count1 < 20
    psp = ps - stepsize*(hh\gg);
    if check_con(psp)==0
        stepsize = stepsize/2; 
        count1 = count1 + 1;
        continue
    end
    dd = norm(gg);
    ps = psp;
    [~,gg,hh] = Grad_Hess_eprior(ps,vSyy,vSxx,vSxy,Eh1sq,T,r);
    count1 = count1 + 1;
end  
[I2new] = Grad_Hess_eprior(ps,vSyy,vSxx,vSxy,Eh1sq,T,r);
phi = ps(1); s = ps(2); 
dQ1 = max(0,sum(I2new) - sum(I2old)); 
end

function ind = check_con(psp)
r = length(psp)/2;
phi = psp(1:r);
s = psp(r+1:end);
if and(sum(abs(phi)<.99)==r,sum(s>0.001)==r)
    ind=1;
else
    ind=0;
end
end 
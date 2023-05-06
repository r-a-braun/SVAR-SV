function [val,grad,hess] = Grad_Hess_eprior(para,vSyy,vSxx,vSxy,Eh1sq,T,r)
phi = para(1:r); 
s = para(r+1:end);
%% Prepare empty matrices  
grad_corr = zeros(2*r,1);
Hess_corr = zeros(2*r);
for i = 1:r
    [grad_ci, Hess_ci] = gradhess_correction( [phi(i);s(i)] , T);
    grad_corr(i:r:end) = grad_ci;
    Hess_corr(i:r:end,i:r:end) = Hess_ci; 
end 

sig2mu = s./(1-phi.^2).*(T*(1-phi.^2)-2*phi.*(1-phi.^T))./(T^2*(1-phi).^2); 

gphi1 = -phi./(1-phi.^2) + phi./s.*Eh1sq;
gphi = 1./(s).*vSxy - 1./s.*vSxx.*phi  + gphi1;


gs1 = - 1./(2.*s) + ((1-phi.^2)./(2.*s.^2)).*Eh1sq;
gs = -(T-1)./2*(1./(s)) + 1./(2.*s.^2).*(vSyy - 2*vSxy.*phi  + phi.^2.*vSxx)  + gs1; 
grad = [gphi;gs]+grad_corr; 


Hphi = (- 1./s.*vSxx ) + (-(1+phi.^2)./((1-phi.^2).^2)  + Eh1sq./s );
Hphis= (-1./(s.^2).*vSxy + 1./(s.^2).*vSxx.*phi) +  (-phi./(s.^2).*Eh1sq);
Hs = ( (T-1)./2*(1./(s.^2)) - 1./(s.^3).*(vSyy - 2*vSxy.*phi  + phi.^2.*vSxx) )...
    + (1./(2.*s.^2)-(1-phi.^2)./(s.^3).*Eh1sq); 
hess = diag([Hphi;Hs])+diag(Hphis,-r)+diag(Hphis,r)+Hess_corr; 
v1 = -0.5*log( s./(1-phi.^2)) -.5*(1./(s./(1-phi.^2))).*Eh1sq;

%val = sum( v1  -0.5.*(T-1).*log(s)  -0.5.*(1./s).*(vSyy - 2.*phi.*vSxy  + phi.^2.*vSxx) );

val = sum( v1  -0.5.*(T-1).*log(s)  -0.5.*(1./s).*(vSyy - 2.*phi.*vSxy  + phi.^2.*vSxx) + .5*log(sig2mu));



end



function [grad,Hess ] = gradhess_correction( phis , T)
%GRADHESSSIG2 Summary of this function goes here
%   Detailed explanation goes here

phi = phis(1);
s = phis(2);
sig2 = s./(1-phi^2)*(T*(1-phi^2)-2*phi*(1-phi^T))/(T^2*(1-phi)^2);
g = zeros(2,1);
H = zeros(2);
id = 1:T-1;

%% Gradients:
ds = 1./(T*(1-phi^2)) + 2./(T^2)*sum( (T-id).*phi.^(id)./(1-phi^2));
dphi = s.*(   2.*phi./(T.*(1-phi.^2).^2) + 2./(T^2) *sum( (T-id).*( id.*phi.^(id-1)./(1-phi.^2) + 2.*phi.^(id+1)./((1-phi.^2).^2))));
ds2 = 0;
dphis =  2.*phi./(T.*(1-phi.^2).^2) + 2./(T^2) *sum( (T-id).*( id.*phi.^(id-1)./(1-phi.^2) + 2.*phi.^(id+1)./((1-phi.^2).^2)));
dA = 2./(T.*(1-phi.^2).^2) + 8/T*phi^2./( (1-phi^2).^3);
dB = ((id-1).*id.*phi.^(id-2))./(1-phi^2) + 2.*id.*phi.^(id)./((1-phi.^2).^2);
dC = 2.*(id+1).*phi.^id./((1-phi.^2).^2) + 8.*(phi.^(id+2))./((1-phi.^2).^3);
dphi2 = s.*(dA + 2./(T^2).*sum((T-id).*(dB+dC)));

g(1) = dphi;
g(2) = ds;
H(1,1) = dphi2;
H(1,2) = dphis;
H(2,1) = dphis;
H(2,2) = ds2;


%% Chain Rule:
dy = .5*1./sig2; 
grad = dy.*g; 
dy2 = -.5*1./(sig2.^2); 
Hess = dy2.*g*g' + dy.*H;
 
end





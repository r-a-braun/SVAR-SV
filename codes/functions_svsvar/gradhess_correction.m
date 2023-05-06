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


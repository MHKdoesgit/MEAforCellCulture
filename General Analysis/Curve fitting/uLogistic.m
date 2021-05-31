function [f,J] = uLogistic(params,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
%unwrap inputs
x=x(:);
m=params(1); 
M1=params(2); x1=params(3); k1=params(4);
M2=params(5); x2=params(6); k2=params(7);
%==========================================================================
%calculate function value
denom1=1+exp((x1-x).*k1);
denom2=1+exp((x2-x).*k2);
f =  m+(M1-m)./denom1+(M2-m)./denom2; 
%==========================================================================
%calculate Jacobian
if nargout>1
    
    J1=1-1./denom1-1./denom2; 
    
    J2=1./denom1;
    J3=-k1*exp((x1-x)*k1).*(M1-m)./denom1.^2;
    J4=(x-x1).*exp((x1-x)*k1).*(M1-m)./denom1.^2;
    
    J5=1./denom2;
    J6=-k2*exp((x2-x)*k2).*(M2-m)./denom2.^2;
    J7=(x-x2).*exp((x2-x)*k2).*(M2-m)./denom2.^2;
    
    J=[J1 J2 J3 J4 J5 J6 J7];
%==========================================================================

end


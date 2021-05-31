function [ F, J ] = logistic4( params, x )
%LOGISTICFUN evaluates modified logistic with params at x
%   Detailed explanation goes here

x=x(:);

F =  params(1)+(params(2)-params(1))./(1+exp((params(3)-x).*params(4))); 

J1=1-1./(1+exp((params(3)-x)*params(4)));
J2=1./(1+exp((params(3)-x)*params(4)));
J3=-params(4)*exp((params(3)-x)*params(4)).*(params(2)-params(1))./(1+exp((params(3)-x)*params(4))).^2;
J4=-(params(3)-x).*exp((params(3)-x)*params(4)).*(params(2)-params(1))./(1+exp((params(3)-x)*params(4))).^2;
J=[J1 J2 J3 J4];


end


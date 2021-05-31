function [params, fcn] = fitULogisticToSpikes(xvals, spikes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
%unwrap and standardize variables
spikes=spikes(:); xvals=xvals(:);
mx=mean(xvals); sx=std(xvals);
xx=(xvals-mx)/sx; yy=spikes/max(spikes);
%==========================================================================
%guess params by fitting single logistic
gparams=marmoset.fitGenLogisticToSpikes(xx,yy);
mguess=gparams(1);
if gparams(4)>0
    onguess=gparams(2:4); offguess=zeros(size(onguess));
    onb=1;
else
    offguess=gparams(2:4); onguess=zeros(size(offguess));
    onb=-1;
end
guess=[mguess onguess' offguess'];
%==========================================================================
%Setting up bounds
lb=[ 0   0  -Inf   0   0  -Inf -Inf]; 
ub=[Inf Inf  Inf  Inf Inf  Inf    0];
%==========================================================================
%Setting up linear inequality constraints
A=[1,-1,0,0,0,0,0;... %m should be smaller than M1
    1,0,0,0,-1,0,0;...%m should be smaller than M2
    0,-onb,0,0,onb,0,0]; %M1/M2 inequality depending on branch
b=[0;0;0]; 
%==========================================================================
%perform nonlinear optimization
foptim=@(p) logisticOptim(p,xx,yy);
[~,startGrad]=foptim(guess); %calculate gradient at initial point
if all(isfinite(startGrad)); useGrad=true; else, useGrad=false; end

options= optimoptions('fmincon','Algorithm','sqp',...
     'Display','off','SpecifyObjectiveGradient',useGrad,...
     'CheckGradients',false);
if isnan(foptim(guess))
    params = NaN(1,7);
    fcn = @(p,x) marmoset.uLogistic(p,x)'; 
    return
end

params = fmincon(foptim, guess,A,b,[],[],lb,ub,[],options);  
fcn = @(p,x) uLogistic(p,x)'; 
%==========================================================================
%bring params back to scale
params([3 6])=params([3 6])*sx+mx; params([4 7])=params([4 7])/sx;
params([1 2 5])=params([1 2 5])*max(spikes);
%==========================================================================
%Xplot = linspace(min(xx), max(xx));
%plot(xx, yy, 'ob',Xplot, fcn(params,Xplot),'-r', Xplot, logistic4(gparams,Xplot),'-k');
end

function [f,g] = logisticOptim(p,X,Y)
%get function value and gradient
[lf,lg] = marmoset.uLogistic(p,X);
%specify optimization problem
f = 0.5*((lf-Y)'*(lf-Y));
if nargout > 1; g = (lf-Y)'*lg; end
end
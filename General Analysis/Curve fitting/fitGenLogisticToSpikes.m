function [params, fcn] = fitGenLogisticToSpikes(xvals, spikes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
%unwrap and standardize variables
spikes=spikes(:); xvals=xvals(:);
mx=mean(xvals); sx=std(xvals);
xx=(xvals-mx)/sx; yy=spikes/max(spikes);
%==========================================================================
yylog=log((yy+1e-3)./(max(yy)-yy+1e-3));
p=polyfit(xx, yylog, 1); %linear fit to get the estimates for the guess

guess=[0; max(yy); -p(2)/p(1); p(1)];

lb=[0 0 -1e3 -1e3];
ub=[1e3 1e3 1e3 1e3];
%==========================================================================
%optimize function
options= optimoptions('fmincon','Algorithm','sqp',...
    'Display','off','SpecifyObjectiveGradient',true);

foptim=@(p) logisticOptim(p,xx,yy);
if isnan(foptim(guess))
    params = NaN(4,1);
    fcn = @(p,x) logistic4(p,x)'; 
    return
end

params = fmincon(foptim, guess,[],[],[],[],lb,ub,[],options);  
fcn = @(p,x) logistic4(p,x)'; 
%==========================================================================
%Xplot = linspace(min(xx), max(xx));
%plot(xx, yy, 'ob',Xplot, logistic4(params,Xplot),'-r',Xplot, logistic4(guess,Xplot),'-k');
%==========================================================================
%rescale params
params(3)=params(3)*sx+mx; params(4)=params(4)/sx;
params([1 2])=params([1 2])*max(spikes);
%==========================================================================
end

function [f,g] = logisticOptim(p,X,Y)
% Calculate objective f

[lf,lg]=logistic4(p,X);

f = 0.5*((lf-Y)'*(lf-Y));

if nargout > 1 % gradient required
    g = (lf-Y)'*lg;
end
    
end
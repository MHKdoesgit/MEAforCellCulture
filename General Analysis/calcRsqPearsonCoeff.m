

function [Rsqr,pearCorrSqrt,explVar] = calcRsqPearsonCoeff(origCurve,predCurve)
%
%%% calcRsqPearsonCoeff %%%
%
%
% This function calcualtes R-square, explained variance and Pearson
% coefficient for the similariy between the original curve and the
% predication curve.
%
%================================Inputs====================================
%
%   origCurve : original data.
%   predCurve : prediction from the model based on the original curve.
%
%================================Output====================================
%
%   Rsqr : R-square for the prediction.
%   pearCorrSqrt : Squared Pearson correlation coefficient for the
%   predicted curve.
%   explVar : Explained vairance for the prediction.
%
% written by Mohammad based on the same function from Helene on 12.06.2018

% checking for correct dimensions
if size(origCurve,1) ~= size(predCurve,1), predCurve = transpose(predCurve); end

sstottemp= origCurve-mean(origCurve);
sstot = sum( (sstottemp.^2));
%regression sum of squares SSreg -> explained sum of squares
ssregtemp= predCurve-mean(origCurve);
ssreg = sum( (ssregtemp.^2));
%sum of squares residuals SSres
ssRestemp=origCurve-predCurve;
ssres = sum( (ssRestemp.^2));
%calculate Rsqr
Rsqr = 1- ssres/sstot;
%calculate expliained variance
explVar= ssreg/sstot;
%corrCoefficient
%https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
corrcoefftemp = corrcoef(origCurve,predCurve);  % cov(X,Y)/std(x)*std(y)

pearCorrSqrt = corrcoefftemp(2).^2;

end
function y = loss_waste_noisefree(theta)
global p
%
% This is the wastewater loss function when noise has been removed (Example 2.8).
%
w=0.1;     	%the weighting coefficient
y = w*(1-theta(1))^2+(1-w)*(1-0.5*theta(1)*theta(2)-theta(2))^2;


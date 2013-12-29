function y=loss_example2_3(theta)
%
% This is the loss function for Styblinski and Tang 
% (1990), example 1.  Used in Example 2.3 in ISSO.  
%
y=.5*(sum(theta.^4)-16*sum(theta.^2)+5*sum(theta));
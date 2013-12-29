function y = example(x)
% Used for numerical integration to determine constant for Example 16.5 on 
% Gibbs sampling with truncated exponential.  Used this function with one of 
% the standard MATLAB numerical integration routines.
y=(ones(1,length(x))-exp(-5*x));
y=y./x
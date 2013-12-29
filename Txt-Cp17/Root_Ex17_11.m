function y = Root_Ex17_10(theta,x)
%This is the derivative of the latest summand in the log-likelihood function for the 
%the organism/dilution problem.
global test truetheta 
% Generate random point z
mu=truetheta;
test=rand;
prob=exp(-mu/x);
if test < prob
   z=0;
else
   z=1;
end
% Compute negative derivative of log-likelihood
y = -(z-1)/x-(z*exp(-theta/x))/(x*(1-exp(-theta/x)));



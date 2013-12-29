% J. C. Spall, January 2002
% Code for doing sequential design in the context of the organism/dilution problem of
% Examples 17.8 and 17.9.  Used in numerical results of Example 17.10.  User needs to 
% provide the gradient input via the "root" function.
clear all
global test truetheta
% meas. noise standard deviation; multplies all elements of N(0,I)noise % vector
a=25;
A=50;
alpha=1; 
n=2000; 		%no. of iterations
cases=20;	%number of cases (replications)
rand('seed',311113)
thetamin=1.59;   %Lower bound on theta 
thetamax=1000;    %Upper bound on theta
root='root_Ex17_10';
%
%the loop 1:cases below is for doing multiple cases for use in %averaging to
%evaluate the relative performance.  
%
errtheta=0;
tolerancetheta=100;	  % Maximum allowable change in an element of theta
truetheta=4;
theta_0=2;
x_0=0.63*theta_0;
for j=1:cases
 caseiter=j
% Initialization 
theta=theta_0;
x=x_0;
%
% The iterations below are the basic R-M iterations. 
%
    for k=0:n-1
       ak=a/(k+1+A)^alpha;
       Yk=feval(root,theta,x);
       thetalag=theta;
       theta=theta-ak*Yk;
       if max(abs(thetalag-theta)) > tolerancetheta;
          theta=thetalag;
       end    
    % Checking for constraints below
       theta=min(theta,thetamax);
       theta=max(theta,thetamin);
       x=0.63*theta;
    end
    errtheta=errtheta+(theta-truetheta)'*(theta-truetheta); 
%  lossthetasq=lossthetasq+loss4thorder(theta)^2;
%  losstheta=losstheta+loss4thorder(theta);
end
% Normalized RMS error over all cases (square root of average MSE) divided by RMS error
% for initial condition
if cases > 1
   ((errtheta/cases)^.5)/(((theta_0-truetheta)'*(theta_0-truetheta))^.5)
else
   errtheta^.5/(((theta_0-truetheta)'*(theta_0-truetheta))^.5)
end
theta
% standard dev. of normalized loss values
%(lossthetasq/(cases*loss4thorder(theta_0)^2)-(losstheta/(cases*loss4thorder(theta_0)))^2)^.5
%losstheta/(cases*loss4thorder(theta_0))

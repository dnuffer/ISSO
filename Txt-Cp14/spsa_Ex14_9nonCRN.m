%J. C. Spall, October 2001
%Some aspects of code tailored to Example 14.9
% on CRN vs. Non-CRN.  This is the non-CRN code.  
global p z lambda;  %declaration of random var. (z) used for normal noise
                 %generation in loss fn. calls given seed above;
                 %also sigma in noise (noise may be dependent on theta)
format long;
p=10;
N=2000;			%total no. of loss measurements (2 x no. of iterations)
cases=50;
alpha =1;
gamma =.166667;
a=.7;
c=.5;          %chosen by standard guidelines
gavg=1;        %no. of grad. estimates averaged in update
A=0;
lambda=[1.1025,1.6945,1.4789,1.9262,0.7505,1.3267,0.8428,0.7247,0.7693,1.3986]';
theta_star=[0.286,0.229,0.248,0.211,0.325,0.263,0.315,0.327,0.323,0.256]';	

lossfinalsq=0;          %variable for cum.(over 'cases')squared loss values
lossfinal=0;            %variable for cum. loss values
theta_lo=zeros(10,1);   %lower bounds on theta  
theta_hi=Inf;           %upper bounds on theta 
tolerancetheta=Inf;     %max. allowable change in any element of theta 
toleranceloss=0;        %tolerance in loss-based blocking step
theta_0=ones(p,1);  
loss='loss_exponential';          % loss function used in algorithm operations
rand('seed',51415927)
losses=zeros(cases,1);%temp statement for use with problem 8.2
thetanorm=zeros(cases,1);%temp statement for use with problem 8.2
for i=1:cases
  theta=theta_0;
  for k=1:N/(2*gavg)
    ak = a/(k+A)^alpha;
    ck = c/k^gamma;
    ghat=0;
    for j=1:gavg
      delta = 2*round(rand(p,1))-1;
      thetaplus = theta + ck*delta;
      thetaminus = theta - ck*delta;
      % Two lines below invoke constraints
      thetaplus=min(thetaplus,theta_hi);
      thetaminus=max(thetaminus,theta_lo);
      yplus=feval(loss,thetaplus);
      yminus=feval(loss,thetaminus);
      ghat = (yplus - yminus)./(2*ck*delta)+ghat;
    end
    thetalag=theta;
    theta=theta-ak*ghat/gavg;
    % Two lines below invoke constraints
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
  end
  lossvalue=0;
  for j=1:p
     lossvalue=lossvalue+lambda(j)/(lambda(j)+theta(j));
  end
  lossvalue=lossvalue+theta'*theta;
  losses(i)=lossvalue; %used in Example 14.9
  lossfinalsq=lossfinalsq+lossvalue^2;
  lossfinal=lossfinal+lossvalue;
  thetanorm(i)=((theta-theta_star)'*(theta-theta_star))^0.5;
  theta  
end
% Display results: Mean loss value and standard deviation
%
disp('mean loss value over "cases" runs') 
lossfinal/cases
%
if cases > 1 
	disp('sample standard deviation of mean loss value') 
   sd=((cases/(cases-1))^.5)*(lossfinalsq/cases-(lossfinal/cases)^2)^.5;
   sd=sd/(cases^.5)
else
end
losses %for use with example 14.9
thetanorm/((theta_0-theta_star)'*(theta_0-theta_star))^0.5 %for use with example 14.9

% J. C. Spall, February 2001
% This code calculates the sample mean of scalar random variables using basic 
% root-finding SA(Robbins-Monro SA).  Also computes confidence interval for sample
% mean
%
n=100000;	%no. of iterations
rand('seed',71111113)
a=4;
A=0;
alpha=1;
theta_0=0;
thetaavg=0;
thetaavgsq=0;
Navg=20;
t=2.093; % value for t-distribution (Navg-1 d.o.f.)
for i=1:Navg
  theta=theta_0;
  for k=1:n
     x=2*rand;
     theta=theta-(a/(k+A)^alpha)*(theta-x);
  end
  thetaavg=(i-1)*thetaavg/i+theta/i;
  thetaavgsq=(i-1)*thetaavgsq/i+(theta^2)/i;
end
mean=thetaavg
stdevmean=(((Navg/(Navg-1))*(thetaavgsq-mean^2))^.5)/(Navg^.5)
confint=[mean-t*stdevmean,mean+t*stdevmean]
   

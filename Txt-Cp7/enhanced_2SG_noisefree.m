% J.C. Spall, Feb. 2006
% enhanced_2SG_noisefree 
% Code for evaluation of enhanced second-order SA, which includes non-identical 
% weighting for the Hessian inputs and feedback.  This code compares the enhanced method 
% with the standard 2SG method in Spall (2000). 
% Code is for comparative evaluation purposes; hence,  
% it includes much that is not required for a basic implementation. This
% code works with noise-free gradient
% measurements, including the simple choice of weighting w_k=w/k^delta. 
% Code allows for checking for simple constraint 
% violation (componentwise constraints). Some of the statements below are specific to the
% 4th-order loss function used in Example 6.6 of ISSO ('loss4thgrad').   
% 
% We use thetaH for the standard 2SG recursion and theta for the enhanced 2SG recursion.   
%
clear all
close all
global p; 
p=10;
B=triu(ones(p,p))/p; %same B as in loss functions for use in true Hessian 
% gain sequence a_k numerator, stability constant, and c_k numerator
a=2;
A=50;
c=.01;
alpha=.602;         %a_k decay rate for 2SPSA
gamma=.101;         %c_k decay rate for 2SPSA
w=1/10;         %numerator in weight sequence for Hessian averaging
d=.501;         %decay rate in weight sequence 
n=100;
loss='loss4thorder'; %loss function for evaluation of algorithm (no noise)
grad='loss4thgrad_noisefree'; %gradient of loss function (no noise) 
% number of cases (replications) of standard and enhanced 2SG methods
cases=50;
tolerancetheta=1;		%max. allowable change in elements of theta
toleranceloss=-100;     %large negative number=no blocking
%avg=1;                  %no. of loss evals. per loss-based blocking step 
rand('seed',31415297)
%randn('seed',311113)
thetamin=-10*ones(p,1);   %Lower bounds on theta (in unconstrained, set to large neg. no.) 
thetamax=10*ones(p,1);    %Upper bounds on theta (in unconstrained, set to large po
%
%the loop 1:cases below is for doing multiple cases for use in averaging to
%evaluate the relative performance.  
%
%
%lines below initialize various recursions for the gradient/Hess. averaging
%and for final error reporting based on the average of the solutions for 
%"cases" replications.
meanHbar=0;
errtheta=0;
errthetaH=0;
losstheta=0;
lossthetaH=0;
lossthetasq=0;
lossthetaHsq=0;
truetheta=zeros(p,1);
theta_0=.2*ones(p,1);
Hhat=eye(p); %dummy statement (sets dimension)
norm=zeros(cases,1); %vector for storing the norms of diff. between 2SG estimated and true Hessian
norm_enh=zeros(cases,1); %vector for storing the norms of diff. between enhanced 2SG estimated and true Hessian
tempH11=zeros(cases,1);
tempH11_enh=zeros(cases,1);
for j=1:cases
  caseiter=j
% Initialization
  theta=theta_0;
  thetaH=theta;
  Hbar=1.05*2*B'*B;     %Not relevant unless setting prior on H for standard 2SG (see also line 95)
  EHbar=1.05*2*B'*B;    %.1*eye(p);%1.05*2*B'*B;      %relevant when w < 1 
  EHbarbar=1.05*2*B'*B; %enhanced Hbarbar i.c.
  %AA=randn(p,p);
  %EHbarbar=2*B'*B %+.0000000001*(AA+AA') %TEMP STATEMENT
  %cumcksq=0; %variable for storing the sum of c_k squared for use in weighting
%
%*********2SG Standard Method********* 
%   
  lossold=feval(loss,theta_0); %lossold for use in loss-based blocking
  for k=0:n-1
    ak=a/(k+1+A)^alpha;
    ck=c/(k+1)^gamma;
    ghatinput=0;
    Hhatinput=0;
% Generation of gradient and Hessian estimates                 
    delta=2*round(rand(p,1))-1;
    thetaplus=thetaH+ck*delta;
    thetaminus=thetaH-ck*delta;  
    ghatplus=feval(grad,thetaplus); 
    ghatminus=feval(grad,thetaminus);
    ghatinput=feval(grad,thetaH); 
    deltaghat=ghatplus-ghatminus;
    for i=1:p
      Hhat(i,:)=deltaghat(i)./(2*ck*delta)';
    end
    Hhatinput=.5*(Hhat+Hhat');
    Hbar=((k+1)/(k+2))*Hbar+Hhatinput/(k+2); %If including prior on H, set k-->k+1        
%   Form below uses Gaussian elimination
    Hbarbar=sqrtm(Hbar*Hbar+.00000001*exp(-k)*eye(p));
    thetaHlag=thetaH;
    thetaH=thetaH-ak*(Hbarbar\ghatinput); 
%   Checking for constraints below
    thetaH=min(thetaH,thetamax); 
    thetaH=max(thetaH,thetamin);    
    lossnew=feval(loss,thetaH);
    if lossnew > lossold-toleranceloss;
      thetaH=thetaHlag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end 
    if max(abs(thetaHlag-thetaH)) > tolerancetheta;
      lossold=lossold1;
      thetaH=thetaHlag;
    end    
 %   thetaH
  end
  %
%*********Enhanced 2SG Method********* 
%  
  lossold=feval(loss,theta_0); %lossold for use in loss-based blocking
  for k=0:n-1
    ak=a/(k+1+A)^alpha;
    ck=c/(k+1)^gamma;
    wk=w/(k+1)^d;
    ghatinput=0;
    Hhatinput=0;
% Generation of gradient and Hessian estimates                 
    delta=2*round(rand(p,1))-1;
    %theta=truetheta; %TEMP STATEMENT
    thetaplus=theta+ck*delta;
    thetaminus=theta-ck*delta;  
    ghatplus=feval(grad,thetaplus); 
    ghatminus=feval(grad,thetaminus);
    ghatinput=feval(grad,theta);
    deltaghat=ghatplus-ghatminus;
    for i=1:p
      Hhat(i,:)=deltaghat(i)./(2*ck*delta)';
    end
    Dk=delta*(1./delta)'-eye(p);
    %EHbarbar=2*B'*B; %TEMP STATEMENT
    Hhatinput=.5*(Hhat+Hhat')-.5*(EHbarbar*Dk+Dk'*EHbarbar);%temp--Hhatinput  
   % cumcksq=cumcksq+ck^2;
    EHbar=(1-wk)*EHbar+wk*Hhatinput; 
    %   Form below uses Gaussian elimination
    EHbarbar=sqrtm(EHbar*EHbar+.00000001*exp(-k)*eye(p));
    %EHbarbar=2*B'*B; %TEMP STATEMENT
    thetalag=theta;
    theta=theta-ak*(EHbarbar\ghatinput);
    % Checking for constraints below
    theta=min(theta,thetamax);
    theta=max(theta,thetamin);
%   Steps below perform "blocking" step with "avg" no. of loss evaluations
    lossnew=feval(loss,theta);
    if lossnew > lossold-toleranceloss;
      theta=thetalag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end 
    if max(abs(thetalag-theta)) > tolerancetheta;
      lossold=lossold1;
      theta=thetalag;
    end
 %   theta
  end
  %loss4thorder(theta)/loss4thorder(theta_0)
  %theta
  %EHbar
  %loss4thorder(theta)
  %meanHbar=meanHbar+Hbar;
  %thetaH
  %loss4thorder(thetaH)/loss4thorder(theta_0) %temp
  errthetaH=errthetaH+(thetaH-truetheta)'*(thetaH-truetheta);
  errtheta=errtheta+(theta-truetheta)'*(theta-truetheta); 
  lossthetaHsq=lossthetaHsq+feval(loss,thetaH)^2;
  lossthetasq=lossthetasq+feval(loss,theta)^2;
  lossthetaH=lossthetaH+feval(loss,thetaH);
  losstheta=losstheta+feval(loss,theta);
  norm(j)=((max(eig((Hbarbar-2*B'*B)*(Hbarbar-2*B'*B))))^0.5)/.895321; %.895321 is norm of H*=2*B'*B
  norm_enh(j)=((max(eig((EHbarbar-2*B'*B)*(EHbarbar-2*B'*B))))^0.5)/.895321;
  %tempH11(j)=Hbar(1,1); %temporary statement to look at variability of H_11 estimate
  %tempH11_enh(j)=EHbar(1,1); %temporary statement to look at variability of H_11 estimate
end
%meanHbar/cases
% normalized results of standard 2SG and enhanced 2SG
norm_thetaH=((errthetaH/cases)^.5)/((theta_0-truetheta)'*(theta_0-truetheta))^.5
norm_theta=((errtheta/cases)^.5)/((theta_0-truetheta)'*(theta_0-truetheta))^.5
% standard dev. of mean of normalized loss values; these are by multiplied by 
% (cases/(cases-1))^.5 to account for loss of degree of freedom in standard 
% deviation calculation before using with t-test
if cases > 1;
    (((cases/(cases-1))^.5)*(lossthetaHsq/(cases*feval(loss,theta_0)^2)-(lossthetaH/(cases*feval(loss,theta_0)))^2)^.5)/(cases^.5)
    (((cases/(cases-1))^.5)*(lossthetasq/(cases*feval(loss,theta_0)^2)-(losstheta/(cases*feval(loss,theta_0)))^2)^.5)/(cases^.5)
end 
norm_lossthetaH=lossthetaH/(cases*feval(loss,theta_0))
norm_losstheta=losstheta/(cases*feval(loss,theta_0))
%Hbar
%EHbarbar

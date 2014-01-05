% J.C. Spall, April 2006
% enhanced_2SPSA_noisefree 
% Code for evaluation of enhanced second-order SA, which includes non-identical 
% weighting for the Hessian inputs and feedback.  This code compares the enhanced method 
% with the standard 2SPSA method in Spall (2000). 
% Code is for comparative evaluation purposes; hence,  
% it includes much that is not required for a basic implementation. This
% code works with noise-free loss
% measurements, including the simple choice of weighting w_k=w/k^delta. 
% Code allows for checking for simple constraint 
% violation (componentwise constraints). Some of the statements below are specific to the
% 4th-order loss function used in Example 6.6 of ISSO ('loss4thorder').  
% 
% We use thetaH for the standard 2SPSA recursion and theta for the enhanced 2SPSA recursion.   
%
clear all
close all
global p;       
p=10;
a=1;           %value of numerator in a_k in 2SPSA part
A=50;           %stability constant 
c=.01;          %numerator in c_k 
ctilda=c;       %numerator in ctilda_k for 2SPSA;
alpha=1;     %a_k decay rate for 2SPSA
gamma=.166667;     %c_k decay rate for 2SPSA
w=1/10;         %numerator in weight sequence for Hessian averaging
d=.501;         %decay rate in weight sequence 
n=5000;        %no. of iterations
loss='loss4thorder';    %loss function for evaluation of algorithm (no noise)
cases=50;               %number of cases (replications)
tolerancetheta=1;		%max. allowable change in elements of theta
toleranceloss=-100;     %large negative number=no blocking
%avg=1;                 %no. of loss evals. per loss-based blocking step 
rand('seed',31415297)
%randn('seed',3111113)
thetamin=-10*ones(p,1);   %Lower bounds on theta (in unconstrained, set to large neg. no.) 
thetamax=10*ones(p,1);    %Upper bounds on theta (in unconstrained, set to large pos. no.)2*B'*B
%sigma=2;                        
%
%lines below initialize various recursions for the gradient/Hess. averaging
%and for final error reporting based on the average of the solutions for 
%"cases" replications.
%
meanHbar=0;
errtheta=0;
errthetaH=0;
lossthetaH=0;       %thetaH corresponds to standard 2SPSA
losstheta=0;        %theta corresponds to enhanced
lossthetaHsq=0;		%cum. sum of loss squared values
lossthetasq=0;
tracelamdaH=0;
tracelamda=0;
truetheta=zeros(p,1);
theta_0=.2*ones(p,1); 
B=triu(ones(p,p))/p;        %matrix used as part of 4th-order loss function
Hbar_0=0*1.05*2*B'*B;       %Not relevant unless setting prior on H for standard 2SG (see also line 110)
Hbarbar_0=0*1.05*2*B'*B;    %Does not enter calculations except with n=0
EHbar_0=1.05*2*B'*B;      %.1*eye(p);%1.05*2*B'*B;      %relevant when w < 1 
EHbarbar_0=1.05*2*B'*B;
Hhat=eye(p); %dummy statement (sets dimension)
norm=zeros(cases,1); %vector for storing the norms of diff. between 2SG estimated and true Hessian
norm_enh=zeros(cases,1); %vector for storing the norms of diff. between enhanced 2SG estimated and true Hessian
%tempH11=zeros(cases,1);
%tempH11_enh=zeros(cases,1);
for j=1:cases
  caseiter=j
% Initialization
  theta=theta_0;
  thetaH=theta;
  Hbar=Hbar_0;     %Not relevant unless setting prior on H for standard 2SG (see also line 103)
  Hbarbar=Hbarbar_0;  %Does not enter calculations except with n=0
  EHbar=EHbar_0;    %.1*eye(p);%1.05*2*B'*B;      %relevant when w < 1 
  EHbarbar=EHbarbar_0; %enhanced Hbarbar i.c.
 % temp=zeros(n,2);
%
% START 2SPSA ITERATIONS FOLLOWING INITIALIZATION
%  
  lossold=loss4thorder(theta_0); %lossold for use in loss-based blocking
  for k=0:n-1
    ak=a/(k+1+A)^alpha;
    ck=c/(k+1)^gamma;
    ctilda_k=ctilda/(k+1)^gamma;
    wk=w/(k+1)^d; 
    ghatinput=0;
    Hhatinput=0;
% Generation of gradient and Hessian estimates 
    delta=2*round(rand(p,1))-1;
    thetaplus=thetaH+ck*delta;
    thetaminus=thetaH-ck*delta;  
    yplus=feval(loss,thetaplus);
    yminus=feval(loss,thetaminus);
    ghat=(yplus-yminus)./(2*ck*delta);  
% Generate perturbed theta values for Hessian update
    deltatilda=2*round(rand(p,1))-1;
    thetaplustilda=thetaplus+ctilda_k*deltatilda;
    thetaminustilda=thetaminus+ctilda_k*deltatilda;
% LOSS FUNCTION CALLS      
    yplustilda=feval(loss,thetaplustilda);
    yminustilda=feval(loss,thetaminustilda);
    ghatplus=(yplustilda-yplus)./(ctilda_k*deltatilda);
    ghatminus=(yminustilda-yminus)./(ctilda_k*deltatilda);
% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS. PER ITERATION
    %ghatinput=((m-1)/m)*ghatinput+ghat/m;
    deltaghat=ghatplus-ghatminus;
    for i=1:p
      Hhat(i,:)=deltaghat(i)./(2*ck*delta)';
    end
    Hhatinput=.5*(Hhat+Hhat');
    %temp(k+1,1)=Hhatinput(1,1); 
    Hbar=(k/(k+1))*Hbar+Hhatinput/(k+1);%No Prior: Hbar=(k/(k+1))*Hbar+Hhatinput/(k+1); if including prior on H, set Hbar=((k+1)/(k+2))*Hbar+Hhatinput/(k+2)
    %Hbar=(1-wk)*Hbar+wk*Hhatinput; %temp
%   THE THETA UPDATE (FORM BELOW USES GAUSSIAN ELIMINATION TO AVOID DIRECT 
%   COMPUTATION OF HESSIAN INVERSE)
    Hbarbar=sqrtm(Hbar*Hbar+.00000001*exp(-k)*eye(p));%Hbarbar(1,1)
    %temp1(k+1,1)=Hbarbar(1,1); %TEMP
    thetaHlag=thetaH;
%   The main update step
    thetaH=thetaH-ak*(Hbarbar\ghat);
%   Two lines below invoke constraints
    thetaH=min(thetaH,thetamax);
    thetaH=max(thetaH,thetamin);
%   Steps below perform "blocking" step with "avg" no. of loss evaluations
    lossnew=feval(loss,thetaH);
    if lossnew > lossold-toleranceloss;
      thetaH=thetaHlag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end 
    if max(abs(thetaHlag-thetaH)) > tolerancetheta; 
      thetaH=thetaHlag;
      lossold=lossold1;
    end 
  end
%  thetaH
%
% ENHANCED 2SPSA
%
  lossold=loss4thorder(theta_0); %lossold for use in loss-based blocking
  for k=0:n-1 
    ak=a/(k+1+A)^alpha;
    ck=c/(k+1)^gamma;
    ctilda_k=ctilda/(k+1)^gamma;
    wk=w/(k+1)^d;
    Hhatinput=0;
% Generation of gradient and Hessian estimates 
    delta=2*round(rand(p,1))-1;
    thetaplus=theta+ck*delta;
    thetaminus=theta-ck*delta;  
    yplus=feval(loss,thetaplus);
    yminus=feval(loss,thetaminus);
    ghat=(yplus-yminus)./(2*ck*delta);  
% Generate perturbed theta values for Hessian update
    deltatilda=2*round(rand(p,1))-1;
    thetaplustilda=thetaplus+ctilda_k*deltatilda;
    thetaminustilda=thetaminus+ctilda_k*deltatilda;
% LOSS FUNCTION CALLS      
    yplustilda=feval(loss,thetaplustilda);
    yminustilda=feval(loss,thetaminustilda);
    ghatplus=(yplustilda-yplus)./(ctilda_k*deltatilda);
    ghatminus=(yminustilda-yminus)./(ctilda_k*deltatilda);
% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS. PER ITERATION
    %ghatinput=((m-1)/m)*ghatinput+ghat/m;
    deltaghat=ghatplus-ghatminus;
    for i=1:p
      Hhat(i,:)=deltaghat(i)./(2*ck*delta)';
    end
    Dk=delta*(1./delta)'-eye(p);
    Dtilda=deltatilda*(1./deltatilda)'-eye(p);
   %EHbarbar=2.01*B'*B; %temp
    Psi_k=Dtilda'*EHbarbar*Dk+Dtilda'*EHbarbar+EHbarbar*Dk;%EHbarbar
    Psi_k=.5*Psi_k+.5*Psi_k'; 
   % Psi_k=0; %temp
    Hhatinput=.5*(Hhat+Hhat')-Psi_k; %Hhatinput(1,1)%, Hhatinput(10,10) 
    %temp(k+1,2)=Hhatinput(1,1); %TEMP 
    %Hhatinput=.5*(Hhat+Hhat'); Hhatinput %temp (reactivate above line)
    EHbar=(1-wk)*EHbar+wk*Hhatinput; %EHbar(1,1)     %temp   
   % EHbar=(k/(k+1))*EHbar+Hhatinput/(k+1);  %temp--turn this off
%   THE THETA UPDATE (FORM BELOW USES GAUSSIAN ELIMINATION TO AVOID DIRECT 
%   COMPUTATION OF HESSIAN INVERSE)
    EHbarbar=sqrtm(EHbar*EHbar+.00000001*exp(-k)*eye(p)); 
%    EHbarbar=EHbar;EHbarbar(1,1)%,EHbarbar(10,10) %temp
    %temp1(k+1,2)=EHbarbar(1,1); %TEMP 
    thetalag=theta;
%   The main theta update step
    theta=theta-ak*(EHbarbar\ghat);
%   Two lines below invoke constraints
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
      theta=thetalag;
      lossold=lossold1;
    end
  end
  %theta
  errthetaH=errthetaH+(thetaH-truetheta)'*(thetaH-truetheta);   %Sum of error in thetaH values
  errtheta=errtheta+(theta-truetheta)'*(theta-truetheta);       %Sum of error in theta values
  lossthetaHsq=lossthetaHsq+feval(loss,thetaH)^2;               %Sum of squared L(thetaH)values
  lossthetasq=lossthetasq+feval(loss,theta)^2;                  %Sum of squared L(theta)values
  lossthetaH=lossthetaH+feval(loss,thetaH);                     %Sum of L(thetaH)values
  losstheta=losstheta+feval(loss,theta);                        %Sum of L(theta)values 
  tracelamdaH=tracelamdaH+trace((Hbarbar-2*B'*B)*(Hbarbar-2*B'*B));%Sum of trace error in "regular" Hessian estimate
  tracelamda=tracelamda+trace((EHbarbar-2*B'*B)*(EHbarbar-2*B'*B));%Sum of trace error in enhanced Hessian estimate
  %B=triu(ones(p,p))/p;
  norm(j)=((max(eig((Hbarbar-2*B'*B)*(Hbarbar-2*B'*B))))^0.5)/.895321; %.895321 is norm of H*=2*B'*B for loss4thorder
  norm_enh(j)=((max(eig((EHbarbar-2*B'*B)*(EHbarbar-2*B'*B))))^0.5)/.895321;
end
% normalized results of standard 2SPSA and enhanced 2SPSA
norm_thetaH=((errthetaH/cases)^.5)/((theta_0-truetheta)'*(theta_0-truetheta))^.5
norm_theta=((errtheta/cases)^.5)/((theta_0-truetheta)'*(theta_0-truetheta))^.5
%if norm(theta_0-truetheta)~= 0
%  ((errtheta/cases)^.5)/norm(theta_0-truetheta)
%end  
% standard dev. of mean of normalized loss values; these are by multiplied by 
% (cases/(cases-1))^.5 to account for loss of degree of freedom in standard 
% deviation calculation before using with t-test
norm_lossthetaH=lossthetaH/(cases*loss4thorder(theta_0));
norm_losstheta=losstheta/(cases*loss4thorder(theta_0));
% Test statistic (t-value) comparing loss values;
% standard dev. of mean of normalized loss values are multiplied by 
% (cases/(cases-1)) to account for loss of degree of freedom in standard 
% deviation calculation before using with t-test
if cases > 1;
    %(((cases/(cases-1))^.5)*(lossthetaHsq/(cases*loss4thorder(theta_0)^2)-(lossthetaH/(cases*loss4thorder(theta_0)))^2)^.5)/(cases^.5)
    %(((cases/(cases-1))^.5)*(lossthetasq/(cases*loss4thorder(theta_0)^2)-(losstheta/(cases*loss4thorder(theta_0)))^2)^.5)/(cases^.5)
    s2_standard=(cases/(cases-1))*(lossthetaHsq/(cases*loss4thorder(theta_0)^2)-(lossthetaH/(cases*loss4thorder(theta_0)))^2);
    s2_enhanced=(cases/(cases-1))*(lossthetasq/(cases*loss4thorder(theta_0)^2)-(losstheta/(cases*loss4thorder(theta_0)))^2);
    ratio_st=s2_standard/cases;
    ratio_en=s2_enhanced/cases;
    stand_dev=(ratio_st+ratio_en)^.5    %this expression is same as (.)^.5 denominator term in (B.4) of ISSO
    dof=-2+(cases+1)*stand_dev^4/((ratio_st)^2+(ratio_en)^2)   %degrees of freedom (top of p. 521 of ISSO)
    t_value=(norm_lossthetaH-norm_losstheta)/stand_dev
end 
tracelamdaH=tracelamdaH/(cases*trace((Hbarbar_0-2*B'*B)*(Hbarbar_0-2*B'*B)))
tracelamda=tracelamda/(cases*trace((EHbarbar_0-2*B'*B)*(EHbarbar_0-2*B'*B)))
norm_lossthetaH
norm_losstheta
% normalized loss values
%losstheta/cases

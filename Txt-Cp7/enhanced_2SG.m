% J.C. Spall, April 2006
% enhanced_2SG
% Code for evaluation of enhanced second-order SA, which includes non-identical 
% weighting for the Hessian inputs and feedback.  This code compares the enhanced method 
% with the standard 2SG method in Spall (2000). 
% Code is for comparative evaluation purposes; hence,  
% it includes much that is not required for a basic implementation. 
% 
% We use thetaH for the standard 2SG recursion and theta for the enhanced 2SG recursion.   
clear all
global z sigma p
p=10;
B=triu(ones(p,p))/p; %same B as in loss functions for use in true Hessian 
% meas. noise standard deviation; multiplies all elements of N(0,I) noise vector
sigma=.05;
% gain sequence a_k numerator, stability constant, and c_k numerator
a=1;
A=0;
c=.05; 
% Exponents in a_k and c_k 
alpha=1; 
gamma=.49;      
n=10000; %no. of iterations
% number of cases (replications) of standard and enhanced 2SG methods
cases=50;
avg=1; % number of noisy loss values to be averaged in blocking step
rand('seed',31415297)
randn('seed',311113)
thetamin=-10*ones(p,1);   %Lower bounds on theta (in unconstrained, set to large neg. no.) 
thetamax=10*ones(p,1);    %Upper bounds on theta (in unconstrained, set to large pos. no.)
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
tolerancetheta=1;
toleranceloss=-100; %large negative number=no blocking
truetheta=zeros(p,1);
theta_0=.2*ones(p,1);
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
  Hbar=0*2*B'*B;     %Not relevant unless setting prior on H for standard 2SG (see also line 90)
  EHbar=0*2*B'*B;    %enhanced Hbar 
  EHbarbar=0*2*B'*B; %enhanced Hbar
  %AA=randn(p,p);
  %EHbarbar=2*B'*B %+.0000000001*(AA+AA') %TEMP STATEMENT
  cumcksq=0; %variable for storing the sum of c_k squared for use in weighting
%
%*********2SG Standard Method********* 
%   
 lossold=0;
  for i=1:avg
    lossold=lossold+loss4thnoise(theta_0);
  end 
  for k=0:n-1
    ak=a/(k+1+A)^alpha;
    ck=c/(k+1)^gamma;
    ghatinput=0;
    Hhatinput=0;
% Generation of gradient and Hessian estimates                 
    delta=2*round(rand(p,1))-1;
    thetaplus=thetaH+ck*delta;
    thetaminus=thetaH-ck*delta;  
    ghatplus=loss4thgrad(thetaplus); 
    ghatminus=loss4thgrad(thetaminus);
    ghatinput=loss4thgrad(thetaH); 
    deltaghat=ghatplus-ghatminus;
    for i=1:p
      Hhat(i,:)=deltaghat(i)./(2*ck*delta)';
    end
    Hhatinput=.5*(Hhat+Hhat');
    Hbar=(k/(k+1))*Hbar+Hhatinput/(k+1); %If including prior on H, set k-->k+1        
%   Form below uses Gaussian elimination
    Hbarbar=sqrtm(Hbar*Hbar+.000001*exp(-k)*eye(p)); 
    thetaHlag=thetaH;
    thetaH=thetaH-ak*(Hbarbar\ghatinput); 
    % Checking for constraints below
    thetaH=min(thetaH,thetamax); 
    thetaH=max(thetaH,thetamin);
%   Steps below perform "blocking" step with "avg" no. of loss evaluations
    lossnew=0;
    for i=1:avg
      lossnew=lossnew+loss4thnoise(thetaH);
    end
    if lossnew/avg > lossold/avg-toleranceloss;
      thetaH=thetaHlag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end 
    if max(abs(thetaHlag-thetaH)) > tolerancetheta;
      lossold=lossold1;
      thetaH=thetaHlag;
    end    
  end
%
%*********Enhanced 2SG Method********* 
%  
  lossold=0;
  for i=1:avg
    lossold=lossold+loss4thnoise(theta_0);
  end
  for k=0:n-1
    ak=a/(k+1+A)^alpha;
    ck=c/(k+1)^gamma;
    ghatinput=0;
    Hhatinput=0;
% Generation of gradient and Hessian estimates                 
    delta=2*round(rand(p,1))-1;
    %theta=truetheta; %TEMP STATEMENT
    thetaplus=theta+ck*delta;
    thetaminus=theta-ck*delta;  
    ghatplus=loss4thgrad(thetaplus); 
    ghatminus=loss4thgrad(thetaminus);
    ghatinput=loss4thgrad(theta);
    deltaghat=ghatplus-ghatminus;
    for i=1:p
      Hhat(i,:)=deltaghat(i)./(2*ck*delta)';
    end
    Dk=delta*(1./delta)'-eye(p);
    Hhatinput=.5*(Hhat+Hhat')-.5*(EHbarbar*Dk+Dk'*EHbarbar); %H estimate with feedback using EHbarbar
    %Hhatinput=.5*(Hhat+Hhat'); %H estimate w/o feedback
    %Hhatinput=.5*(Hhat+Hhat')-.5*(EHbar*Dk+Dk'*EHbar); %H estimate with feedback using EHbar
    cumcksq=cumcksq+ck^2;
    EHbar=((cumcksq-ck^2)/cumcksq)*EHbar+(ck^2)*Hhatinput/cumcksq; 
    EHbarbar=sqrtm(EHbar*EHbar+.000001*exp(-k)*eye(p)); 
    thetalag=theta;
    theta=theta-ak*(EHbarbar\ghatinput);  % Theta update uses Gaussian elimination
    % Checking for constraints below
    theta=min(theta,thetamax);
    theta=max(theta,thetamin);
%   Steps below perform "blocking" step with "avg" no. of loss evaluations
    lossnew=0;
    for i=1:avg
      lossnew=lossnew+loss4thnoise(theta);
    end
    if lossnew/avg > lossold/avg-toleranceloss;
      theta=thetalag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end 
    if max(abs(thetalag-theta)) > tolerancetheta;
      lossold=lossold1;
      theta=thetalag;
    end    
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
  lossthetaHsq=lossthetaHsq+loss4thorder(thetaH)^2;%loss4thorder(thetaH)/loss4thorder(theta_0) %temp
  lossthetasq=lossthetasq+loss4thorder(theta)^2;%loss4thorder(theta)/loss4thorder(theta_0) %temp
  lossthetaH=lossthetaH+loss4thorder(thetaH);
  losstheta=losstheta+loss4thorder(theta);
  norm(j)=(max(eig((Hbarbar-2*B'*B)*(Hbarbar-2*B'*B))))^0.5;
  norm_enh(j)=(max(eig((EHbarbar-2*B'*B)*(EHbarbar-2*B'*B))))^0.5;
 % tempH11(j)=Hbar(1,1); %temporary statement to look at variability of H_11 estimate
 % tempH11_enh(j)=EHbar(1,1); %temporary statement to look at variability of H_11 estimate
end
%meanHbar/cases
% normalized results of standard 2SG and enhanced 2SG
norm_thetaH=((errthetaH/cases)^.5)/((theta_0-truetheta)'*(theta_0-truetheta))^.5
norm_theta=((errtheta/cases)^.5)/((theta_0-truetheta)'*(theta_0-truetheta))^.5
% Test statistic (t-value) comparing loss values;
% standard dev. of mean of normalized loss values are multiplied by 
% (cases/(cases-1)) to account for loss of degree of freedom in standard 
% deviation calculation before using with t-test
norm_lossthetaH=lossthetaH/(cases*loss4thorder(theta_0));
norm_losstheta=losstheta/(cases*loss4thorder(theta_0));
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
norm_lossthetaH
norm_losstheta

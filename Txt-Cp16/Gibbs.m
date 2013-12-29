%J. C. Spall, June 2001
%Performs Gibbs sampling.  Code is set up for Example 16.5 in ISSO (on bounded exponential
%sampling).  Basic logic can be used for other cases.
clear all
B=5;			%Bound for exponential variables
rand('seed',31415927)
cases=1800;
burn=10; 	%No. of iterations for "burn-in" at each replication
post_burn=30;   %No. of iterations for use in computing ergodic averages per replication
erg_mean=zeros(cases,1);	%vector storing the ergodic mean of X(1), r.v. of interest, for each rep.
output=zeros(cases,1);	%vector storing the terminal value for X(1), the r.v. of interest
for i=1:cases
   X(2)=2.5;	%Initial state
   ergodic=0;
   for k=1:burn
      c=1-exp(-B*X(2));
      X(1)=-log(1-c*rand)/X(2);
      c=1-exp(-B*X(1));
      X(2)=-log(1-c*rand)/X(1);
   end
   for k=1:post_burn
      c=1-exp(-B*X(2));
      X(1)=-log(1-c*rand)/X(2);
      c=1-exp(-B*X(1));
      X(2)=-log(1-c*rand)/X(1);
      ergodic=ergodic+X(1);
   end
   erg_mean(i)=ergodic/post_burn;
   output(i)=X(1);
end
%erg_mean
%mean(erg_mean)
%std(erg_mean)
hist(output,20);
       
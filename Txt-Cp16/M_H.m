%J. C. Spall, June 2001
%Performs the M-H sampling.  Code is specialized for Example 16.1 on generating bivariate
%normal r.v.s with mean zero.  Code is not optimized for efficiency (e.g., includes matrix
%inverse).
clear all;
sigma=[1 .9;.9 1];			%covariance matrix
rand('seed',141415927)
cases=3;
burn=500; 	%No. of iterations for "burn-in" at each replication
post_burn=15000;   %No. of iterations for use in computing ergodic averages per replication
%erg_mean=zeros(cases,1);	%vector storing the ergodic mean of X(1), r.v. of interest, for each rep.
%output=zeros(cases,1);	%vector storing the terminal value for X(1), the r.v. of interest
percent_rej=zeros(cases,1); %variable storing the % of rejects for each run
for i=1:cases
   reject=0;
   X=[-1,1]';	%Initial state
   ergodic=0;
   for k=1:burn
      W=X+rand(2,1)-0.5*ones(2,1);	%candidate: addition of U(-.5,.5) r.v.s to current X
      test=rand;
      rho=exp(-.5*W'*inv(sigma)*W)/exp(-.5*X'*inv(sigma)*X);
      if test < rho
         X=W;
      else
         reject=reject+1;
      end
   end
   for k=1:post_burn
      W=X+rand(2,1)-0.5*ones(2,1);	%candidate: addition of uniform r.v.s to current X
      test=rand;
      rho=exp(-.5*W'*inv(sigma)*W)/exp(-.5*X'*inv(sigma)*X);
      if test < rho
         X=W;
      else 
         reject=reject+1;
      end
      ergodic=(k-1)*ergodic/k+(X(1)+X(2))/k;
      erg_plot(k,i)=ergodic;
   end
   percent_rej(i)=100*reject/(burn+post_burn);
end
plot(erg_plot);
percent_rej
%erg_mean
%mean(erg_mean)
%std(erg_mean)
%hist(output,20);
       
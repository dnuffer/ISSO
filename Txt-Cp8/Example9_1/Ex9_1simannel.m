function [THETA, T0, Loss] = simannel()
%by Fred Dilley, spring 2000
%Simulated annealing code.  Uses geometric decay of temperature.
%Provides two ways of dealing with noisy
%loss measurements: one is by using the rho factor to alter the 
%decision criterion and the other is by simple averaging of the loss
%measurements.

randn('seed',1111113)
rand('seed',31415927)

global p 

n=8000;		   %total no. of loss measurements(iterations/lossavg) 
niter=40;		%no. of iters. per temp. setting 	
bk=1;		   	%"Boltzmann's constant"
lambda=.95;		%cooling rate (<=1)

cases = 10;
for i = 1:cases
   
   theta_0=PathPerm(p);
	PathCost(theta_0);

	T=70;			%initial temperature
	theta=theta_0;
	E_old=PathCost(theta);

	for j=1:n/(niter)
   	for k=1:niter
      	perturbTheta = PerturbPath(theta);
	      E_new=PathCost(perturbTheta);
   	   if E_new < E_old
      	   theta=perturbTheta;
         	E_old=E_new;
	      else 
   	      prob=exp(-(E_new-E_old)/(bk*T));
      	   test=rand;
      		if test < prob
         		theta=perturbTheta;
	         	E_old=E_new;
   	   	else
      		end
	   	end
		end
		T=lambda*T;
	end
   
   THETA(i,:) = theta;
   T0(i,:) = theta_0;
   Loss(i,:) = PathCost(theta);
end
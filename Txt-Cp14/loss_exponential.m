function y = loss_exponential(theta)
%This is the loss function used in Example 14.9 (same as Kleinman, et al. 1999). 
global p z lambda
cum=0;
for i=1:p
   z=rand;
   z=-log(1-z)/lambda(i);
   cum=cum+exp(-z*theta(i));
end
y=theta'*theta+cum;
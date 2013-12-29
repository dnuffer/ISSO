function y = loss_exponential_CRN(thetaplus,thetaminus)
%This is the loss function used in Example 14.9 (same as Kleinman, et al. 1999).  This is the CRN version. 
global p z lambda
cum=0;
for i=1:p
   z=rand;
   z=-log(1-z)/lambda(i);
   cum=cum+exp(-z*thetaplus(i))-exp(-z*thetaminus(i));
end
y=thetaplus'*thetaplus-thetaminus'*thetaminus+cum;
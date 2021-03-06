function y=fitpop(thetabit,bits,loss)
% This function returns a vector of N fitness values.
global p N M thetamin thetamax
theta=zeros(p,1);
for i=1:N
   pointer=1;
   for j=1:p
      theta(j)=bit2num(thetabit(i,pointer:pointer+bits(j)-1),...
         bits(j),M(j),thetamin(j),thetamax(j));
      pointer=pointer+bits(j);
   end
   fitness(i)=-feval(loss,theta);
end
y=fitness;
   
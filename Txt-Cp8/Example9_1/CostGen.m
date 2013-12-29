function c = CostGen()
%by Fred Dilley, spring 2000
global p cost
%rand('seed',31415927); %for 10
rand('seed',45939219); %for 100

for i = 1:p
   for j = i+1:p
      cost(i,j) = round(115*rand(1))+10;
      cost(j,i) = cost(i,j);
   end
end
c=cost;
% This program generates the transition matrix for the Modified Genetic Algorithm
clear;
mu = 0.05;  		% this is the mutation probability
chi = 1.0;  		% this is the crossover probability

Npop = 4;   		% this is the population size
lngth = 4;      	% this is the length of the bit string representing the integers
sols = 2^lngth; 	% this is the number of possible bit solutions of lngth 4  

	%this is the number of possible populations
N = factorial(Npop+sols-1)/(factorial(sols-1)*factorial(Npop));  
N
for i = 1:sols
   solution(i,1) = i-1;
end

%solbits = zeros(sols,lngth);

for i = 1:sols
   for j=1:lngth
      k=lngth-j+1;
      solbits(i,k) = bitget(solution(i,1),j);	%solbits has the bit representation of the  
   end														%candidate integer solution
end

%fxsol is the value of the function evaluated at the candidate solutions
[fxsol] = fxnew(solution);

%determine the priority of the solutions 
%fxsolsrt is the function output values sorted in ascending order
%prtysol is a pointer to fxsol telling which position in fxsol is the sorted value in fxsolsrt

[fxsolsrt,prtysol] = sort(fxsol);
k=0;

%go through them in reverse order since we want descending order 
for i = sols:-1:1
   k=k+1;
   %put the actual priority value from 1 to sols in column 1
   priority(k,1) = k;
   j = prtysol(i);
   %put the solution value in the second column of priority 2
   priority(k,2) = solution(j);
end

      
%generate now the number of possible populations for this problem.
%currently I don't have a real general way to do that but here is one that works for a fixed case

%pops is the matrix holding all possible populations of size Npop 3 here 
pops = zeros(Npop,N);

count = 0;

for i = 1:sols   
   for j = i:sols      
      for k = j:sols
         for l = k:sols
         
         count = count+1;
         pops(1,count) = solution(i,1);
         pops(2,count) = solution(j,1);
         pops(3,count) = solution(k,1);
	 		pops(4,count) = solution(l,1);
        
        end
      
    end
   end
end

s = zeros(sols,N);
%calculates the number of each possible solution in each possible population
for i = 1:N
   
   for j = 1:Npop
      
      for k = 1:sols
         
         if  pops(j,i) == solution(k,1) 
             s(k,i) = s(k,i)+1;
          end
          %this part of the loop calculates the priority for each possible solution in each 
          %possible population
          %Note that we only need to save the minimum priority in each population 
          if pops(j,i) == priority(k,2)
             popspty(j,i) = priority(k,1);
          end
      end
   end
end

rstar1 = min(popspty);

%rstar contains the minimum priority value in each population
%pri_rstar contains the sorted minimum priority for each population
%pop_ordr tells which population each came from  
[pri_rstar,pop_ordr] = sort(rstar1);

 for i = 1:N
   rstar(i) = priority(rstar1(i),2); 
   pri_rstar(1,i) = pri_rstar(1,i);
   pri_rstar(2,i) = i;
   pri_rstar(3,i) = pop_ordr(1,i);
end

for i = 1:N
   for j=1:N
      if i==pop_ordr(1,j)
         pristar(i) = j;
      end
   end
   new_s(:,i) = s(:,pop_ordr(1,i));
   new_pops(:,i) = pops(:,pop_ordr(1,i));
   new_rstar(1,i) = rstar(1,pop_ordr(1,i));
end

clear s, pops, rstar; 

R=zeros(sols,sols);
for i =0:sols-1
   for h = 0:sols-1
       i
      ii = i+1;
      hh = h+1;
      R(ii,hh) = RIH(mu, chi, lngth, i, h);
   end;
end;

%function phi = cmpphi(fxsol,s,i,k,sols);
%Need to check this code to see if the indexing is correct


for k = 1:N
   sum=0;
   k
   for jj = 0:sols-1,
      h=jj+1;
      sum = fxsol(h)*new_s(h,k) + sum;
   end
   sstar(1,k) = sum;
	for jj = 0:sols-1
   	  h=jj+1;
      %sum = fxsol(h)*s(h,k) + sum;
      %Y(h,k) = cmpyx(new_s,jj,k,new_rstar,solution);
      Pivalue(h,k) = cmppi(mu,chi,new_s,sstar,fxsol,sols,lngth,R,jj,k);
   end
end





Q=zeros(N,N);
fact = factorial(Npop-1);
for k = 1:N
   
   for v = 1:N
       k
       v
      if pri_rstar(1,v)<= pri_rstar(1,k)      
          prod = 1; 
          for j = 0:sols-1
             jj=j+1;
              Y = cmpyx(new_s,j,v,k,new_rstar,solution);
            % YF = 1/factorial(Y);
            % pivalue = cmppi(mu,chi,s,sstar,fxsol,sols,lngth,R,j,k); 
             prod = prod*(Pivalue(jj,k)^(Y)/factorial(Y)); 
          end
         
          Q(k,v) = fact*prod;
         
       end
       if pri_rstar(1,k)< pri_rstar(1,v)
          Q(k,v) = 0;
       end
       
   end
end
save newrslt_052703_4_4.mat Q N Npop solution new_pops fxsol ;
load newrslt_052703_4_4.mat;
clear sum;
QT = Q';
qx = sum(QT);

q0 = ones(1,N)/N;

%for j = 1:N,
%    for i = 1:N,
%        Q_new(i,j)  = QT(i,j) / qx(j);
%    end;
%end;


%[p,count] = find_p(Q,16);
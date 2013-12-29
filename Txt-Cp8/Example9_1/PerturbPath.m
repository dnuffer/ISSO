function PP = PerturbPath(theta)
%by Fred Dilley, spring 2000
global p

i = round((p-1)*rand(1))+1;
j = round((p-1)*rand(1))+1;

while i == j
 	j = round((p-1)*rand(1))+1;
end

if i > j
  	temp = i;
   i = j;
	j = temp;
end

if rand < .75 %inversion
   
	for m = 1:p
   	if m <= i
      	PP(m) = theta(m);
	   elseif m > j
   	   PP(m) = theta(m);
	   else
   	   PP(m) = theta(i+j+1-m);
	   end
	end
else
   if rand < .5 %Translation 
      if j == p
         PP(1:j-i+1) = theta(i:j);
         PP(j-i+2:p) = theta(1:i-1);
      else
	      for m = 1:p
   	      if m < i
      	      PP(m) = theta(m);
         	elseif m == i
	            PP(m) = theta(j+1);
   	      elseif m > j+1
      	      PP(m) = theta(m);
	         else
   	         PP(m) = theta(m-1);
            end
         end
      end
   else % Switching 
      PP = theta;
      tmp = theta(i);
      PP(i) = PP(j);
      PP(j) = tmp;
   end   
end

      
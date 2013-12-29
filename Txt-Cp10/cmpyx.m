function Y = cmpyx(new_s,j,v,k,new_rstar,solution);
%Add one here to account for calling program starting at 0 and going to sols-1 
j=j+1;

if solution(j) == new_rstar(k)
   if new_s(j,v) > 0 
      Y = new_s(j,v)-1;
      
   else 
   %   Y=0;
      Y=1;
   end
else 
   Y = new_s(j,v); %correct values from Suzuki
end

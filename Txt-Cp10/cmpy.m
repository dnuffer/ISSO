function Y = cmpy(new_s,r,k,new_rstar,solution);
%Add one here to account for calling program starting at 0 and going to sols-1 
r=r+1;
if solution(r) == new_rstar(k)  
 
   Y = new_s(r,k)-1; %correct values from Suzuki
    
else 

   Y = new_s(r,k); %correct values from Suzuki
end

   
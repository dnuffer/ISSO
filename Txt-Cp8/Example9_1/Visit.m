function res = Visit(k)
global id val p cost path paths
%by Fred Dilley, spring 2000
id = id+1;
val(k) = id;

if id == p 
   if path == 0
      path = path + 1;
      paths(1, 1:p) = val';
      paths(1, p+1) = PathCost(val);
   else
      tempCost = PathCost(val);
      if tempCost < paths(1,p+1)
	      paths(1, 1:p) = val';
   	   paths(1, p+1) = PathCost(val);
      end
   end
end

for t = 1:p
   if val(t)==0 
      Visit(t)
   end
end

id = id-1;
val(k)=0;


      
   

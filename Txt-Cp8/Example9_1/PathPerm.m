function PP = PathPerm(p)
%by Fred Dilley, spring 2000
C = zeros(p,1);

for i = 1:p
   C(i) = i;
end

for i = 1:p
   idex = round((p-i)*rand(1))+1;
   PP(i) = C(idex);
   T = zeros(p-i,1);
   k=0;
   for j = 1:p-i+1
      if idex ~= j
         k=k+1;
         T(k) = C(j);
      end
   end
   C = T;
end

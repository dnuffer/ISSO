function Res = RmvExtra(str, Nset)
%by Fred Dilley, spring 2000
global p

Res = zeros(length(str)-length(Nset), 1);

m = 1;
for i = 1:p
   inSet = 0;
   for j = 1:length(Nset)
      if str(i) == Nset(j)
         inSet = 1;
      end
   end
   if inSet == 0
   	Res(m) = str(i);
      m=m+1;
   end
end

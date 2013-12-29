function c = PathCost(l)
%by Fred Dilley, spring 2000
global cost p

c = 0;

for i = 1:p-1
   c = c + cost(l(i),l(i+1));
end

c = c + cost(l(p),l(1));

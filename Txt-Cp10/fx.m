function [y]=fx(x);
%input vector x
jj = size(x,1);

for i = 1:jj
	y(i,1) = -(x(i)-8.0)^2+100;
end;

 

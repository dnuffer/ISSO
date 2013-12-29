function [y]=fx23(x);
%input vector x
jj = size(x,1);

for i = 1:jj
	y(i,1) = abs(x(i)*sin(sqrt(abs(x(i)))))+1;
end;

 

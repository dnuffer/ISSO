function [y]=fxnew(x);
%input vector x
jj = size(x,1);
%x = x-8;
for i = 1:jj
	y(i,1) = -(0.5*((x(i)-4)^4.0-16*(x(i)-4)^2.0+5*(x(i)-4)))+7000000;
end;

 

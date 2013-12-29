function y=bit2num(th_bit,b,m,thetmin,thetmax)
% This function converts bit representation to floating point no. 
% with m decimal places of accuracy (m here represents the relevent
% scalar element of the M vector used in GAbits and fitpop).  This
% applies to one element of the theta vector.
%
unused=2^b-1-(10^m)*(thetmax-thetmin); %this should be 0; included to 
                                       %allow for future handling of 
                                       %general constraints                           
y=0;
for i=1:b
   y=y+th_bit(b+1-i)*2^(i-1);
end
y=thetmin+(thetmax-thetmin)*y/(2^b-1-unused);
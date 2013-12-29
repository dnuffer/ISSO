%J. C. Spall, October 2001
%Computes variance s^2 for CRN version of indifference zone method according to (14.19) in
%ISSO.  Set up below for use in Example 14.11 in ISSO.
%
%no. of possible systems (K) and data points per system (n)
K=4;
n=15;
%data below (cut/paste directly from command window) 
y=[0.6668    0.5700    0.3782    1.6261
    1.4567    0.0196    0.4622   -0.0320
    0.3773   -1.5258    0.8463    0.4510
    0.4182    1.3412    0.2890    0.5391
   -1.6319   -1.6598   -0.9810   -0.1349
   -0.7135   -1.0238    0.6777   -0.4344
    0.0215    0.2789    0.1031    0.7099
   -0.7797   -1.5196   -1.0706   -0.0211
   -1.7013   -0.8083   -1.2284   -1.2433
   -1.0204   -0.3323    0.0745    0.0196
    0.6019   -0.0741   -0.9736    0.2563
    0.5237    0.9471    1.0767    1.7770
   -0.5797    0.2420    0.0483   -0.2027
    1.9999    0.5571   -0.2574    2.1409
    0.4698   -0.3017    0.9283    1.5375]; % data shown are used in Example 14.11 
 %
 %computation of three required mean values below
 M_1=mean(y);
 M_2=mean(y');
 M_3=mean(M_1);
 %double loop to calculate numerator of s^2
 num=0;
 for i=1:K
    for k=1:n
       num=num+(y(k,i)-M_1(i)-M_2(k)+M_3)^2;
    end
 end
 num=2*num;
 variance=num/((K-1)*(n-1))
 
 
 
 
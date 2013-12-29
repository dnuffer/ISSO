function y = loss_waste_noise(theta)
global p noise
sigma1=1;
sigma2=1;
noise(1)=sigma1*randn(1,1);
noise(2)=sigma2*randn(1,1);
w=0.1;     	%the weighting coefficient
y = w*(1-theta(1)+noise(1))^2+(1-w)*(1-0.5*theta(1)*theta(2)-theta(2)+noise(2))^2+1;


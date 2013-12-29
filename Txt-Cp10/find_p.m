function [p,count] =  find_p(Q,k);

N = size(Q,1);

p0 = ones(1,N);
p0 = p0/N;

p=0;
count = 0;
p_temp = p0
QT= Q';

Q5 = QT*QT*QT*QT*QT;
p = p0*Q5

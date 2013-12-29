function R = RIH(mu, chi, lngth, i, h);

%mu  is the mutation rate for the GA
%chi is the crossover rate for the GA
%lngth is the length of the bit string used for the solutions 
%i is just an index representing the ith candidate solution
%h is just an index representing the hth candidate solution
%we will need the bit equivalent of i and h for some computations
%solbits is the bit equivalent of the integers

part1a = (1-mu)^((lngth-sum(bitrep(i,lngth))))*(mu)^(sum(bitrep(i,lngth)));
part1b = (1-mu)^((lngth-sum(bitrep(h,lngth))))*(mu)^(sum(bitrep(h,lngth)));
 
part1 = (part1a+part1b)*0.5*(1-chi);

sum1 = 0;

for nu = 1:lngth-1
   twonum1 = 2^nu-1;
   frst = mu^(sum(bitrep(bitand(twonum1,h),lngth)) + sum(bitrep(i,lngth)) - sum(bitrep(bitand(twonum1,i),lngth)));
   scnd = (1-mu)^(lngth-sum(bitrep(i,lngth)) + sum(bitrep(bitand(twonum1,i),lngth))-sum(bitrep(bitand(twonum1,h),lngth)));
   sum1 = frst * scnd + sum1;
end;
part2 = (chi/(lngth-1))*(sum1);
%part2
sum2 = 0;
for nu = 1:lngth-1
   twonum1 = 2^nu-1;
   frst = mu^(sum(bitrep(bitand(twonum1,i),lngth))+sum(bitrep(h,lngth))-sum(bitrep(bitand(twonum1,h),lngth)));
   scnd = (1-mu)^(lngth - sum(bitrep(h,lngth))+sum(bitrep(bitand(twonum1,h),lngth))-sum(bitrep(bitand(twonum1,i),lngth)));
   sum2 = frst * scnd +sum2;
end
part3 = (chi/(lngth-1))*(sum2);
%part3
 R = part1 + 0.5*(part2+part3);  
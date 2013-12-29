function phi = cmpphi(fxsol,s,sstar,i,k,sols);
%Need to check this code to see if the indexing is correct

%sum=0;
ii=i+1;
%for jj = 0:sols-1
%   h=jj+1;
%   sum = fxsol(h)*s(h,k) + sum;
%end


phi = fxsol(ii)*s(ii,k)/sstar(1,k);



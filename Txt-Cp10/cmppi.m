function pivect = cmppi(mu,chi,s,sstar,fxsol,sols,lngth,R,j,k);
sum=0;
%Need to check on the indexing here
for i = 0:sols-1
    
   for h = 0:sols-1
      ii = i+1;
      hh = h+1;
      iorj = bitxor(i,j);
      horj = bitxor(h,j);
      phi1 = cmpphi(fxsol,s,sstar,iorj,k,sols);
      phi2 = cmpphi(fxsol,s,sstar,horj,k,sols);
      sum = sum + phi1*phi2*R(ii,hh);
      
   end
end
pivect=sum;
if pivect >=1 
   pivect
   j
   k
end

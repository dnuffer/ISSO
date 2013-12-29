function retbits = bitrep(a,lngth);

%This program returns the bit representation of integer a of length lngth stored in column form

for j=1:lngth
      k=lngth-j+1;
      retbits(k,1) = bitget(a,j);%solbits has the bit representation of the  
end														%candidate integer solution

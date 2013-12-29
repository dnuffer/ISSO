function q = compq(Q,N);
%Q=sparse(Q);
D=zeros(N,N);
DM = zeros(N,N);

  for j = 1:N
     j
   k1=j-1;
   k2=j+1;
   NN = N-1;
   if j==1
      D = Q(2:N,2:N);
   elseif j==N
      D = Q(1:NN,1:NN);
   else 
      DM = [Q(1:k1,:);Q(k2:N,:)];
      D = [DM(:,1:k1) DM(:,k2:N)];
   end
   
   Dvect(j) = det(D);
   if Dvect(j) < 0 
      Dvect(j) = 0;
   end
   
end

q = Dvect / sum(Dvect);
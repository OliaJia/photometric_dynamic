A=rand(7,100);
b=rand(7,1);
kmax=5;
err=0.01;
A(:,8)=0;
%[c,k] = Orthogonal_Matching_Pursuit(A, b, kmax, err);
[x_T,indx_set,x,r,normR,residHist, errHist] = OMP( A, b, kmax);
%[cs,i,k] = Orthogonal_Matching_Pursuit(A, b, kmax, err);

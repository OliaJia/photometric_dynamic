function [index]=ReSelectLight(X, NL,last)
[U,S,V]=svd(X,0);

US=U(:,1:3)*S(1:3,1:3);
k=NL;
[index]=greedy_omp_matrix2(X,US,k,last);


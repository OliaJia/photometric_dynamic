function [index]=LightSel(X, NL)
[U,S,V]=svd(X,0);

US=U(:,1:NL)*S(1:NL,1:NL)';
k=NL;
[index]=greedy_omp_matrix(X,US,k);


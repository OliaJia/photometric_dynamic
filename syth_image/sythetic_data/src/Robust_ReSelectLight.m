function [index,E]=Robust_ReSelectLight(X, NL,last)

addpath('E:\Spring2016_Harvard\Research299r\syth_image\fastRPCA-master\fastRPCA-master')
lambda=0.4;
opt.sum=true;
[X,E,~,~] = solver_RPCA_SPGL1(X,lambda,0.001,[],opt);

[U,S,V]=svd(X,0);

UV=U(:,1:3)*S(1:3,1:3)';
k=NL;
[index]=greedy_omp_matrix2(X,UV,k,last);


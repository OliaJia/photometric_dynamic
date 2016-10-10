function [index]=Robust_LightSel(X, NL)

addpath('E:\Spring2016_Harvard\Research299r\syth_image\fastRPCA-master\fastRPCA-master')
lambda=0.4;
opt.sum=true;
[X,E,~,~] = solver_RPCA_SPGL1(X,lambda,0.001,[],opt);
[U,S,~]=svd(X,0);

US=U(:,1:NL)*S(1:NL,1:NL);
k=NL;
[index]=greedy_omp_matrix(X,US,k);


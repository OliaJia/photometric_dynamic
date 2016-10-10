% This program do the column selection
% with the greedy search algorithm 

%% Generate test data
light=7;
pix=2000;


i=randi([1,2000],100,1);
j=randi([1,7]);
v=randi([-5,5],100,1);
Eref=sparse(i,j,v,2000,7);
P1=rand(pix,4);
P2=rand(4,light);
Xr=P1*P2;
Xr=Xr+Eref;
k=3;
addpath('E:\Spring2016_Harvard\Research299r\syth_image\fastRPCA-master\fastRPCA-master')
lambda=0.4;
opt.sum=true;
 [X,E,~,~] = solver_RPCA_SPGL1(Xr,lambda,0.001,[],opt);

[U,S,V]=svd(X,0);

UV=U(:,1:3)*V(:,1:3)';

[ind]=greedy_omp_matrix(X,UV,k);

C=X(:,ind);
error=norm(X-C*inv(C'*C)*C'*X,'fro');

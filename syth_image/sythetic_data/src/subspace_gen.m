%% Generate data to for 2-cluster subspace clustering

function [X,label,L,Norm]=subspace_gen(N,dim,nos_level,sp_level,dist)
% Input:
% N: number of points in space
% dim: dimension of subspaces
% er_level: tge level of error
% distribution of two subspaces

% Output:
% X: the generated dataset

L=rand([dim+1,dim]);
Norm=rand([dim,N]);

L=normc(L')';
Norm=normc(Norm);

X_clean=L*Norm;

index_z=randperm(N);

X_clean(dim+1,index_z(1:floor(N*dist)))=0;

sub1=index_z(1:floor(N*dist));

error=nos_level*rand(dim+1,N);

sp_error=sp_level*sprand(dim+1,N,0.05);

X=X_clean+error+sp_error;

X=abs(X);

label=zeros(N,1);

label(sub1,1)=1;
function [ind,B]=greedy_omp_matrix(X,UV,k)

% Input matrix X to select column in
% k the number of column will be selected
% UV the result of svd
% C the matrix selected, 
% ind the index of column be selected

%% normalized X
N=size(X,1);
M=size(X,2);
A=X;
A=normc(X);
Ap=A;
B=UV;
%res=zeros(k,1);
i=1;
ind=zeros(1,1);
% OMP for matrix search
while i<=k
    prd=B'*Ap;
    [~,index]=max(sum(sqrt(prd.*prd),1));
    B=B-A(:,index)*(B'*A(:,index))';
    ind(i,1)=index;
    A=A-A(:,index)*(A'*A(:,index))';
    A=normc(A);
    Ap=A;
    Ap(:,ind)=0;
    %res(i,1)=norm(Ap,'fro');
    i=i+1;
end



    
    

function [C,Id,missrate,group_real]=SparseSubspace(X,kmax,err,s,rho,w)

N=size(X,2);
% C=zeros(kmax,N);
% Id=zeros(kmax,N);
C=zeros(N,N);
Id=zeros(N,N);
%k=zeros(N,1);
for i=1:1:N
    Xi=X;Xi(:,i)=0;
    %[C(:,i),Id(:,i),k(i,1)]=Orthogonal_Matching_Pursuit(Xi, X(:,i), kmax, err);
    %[C(:,i),Id(:,i),~,~,~,~,~] = OMP(Xi, X(:,i), kmax);
    [C(:,i),Id(:,i),~] = Orthogonal_Matching_Pursuit(Xi, X(:,i), kmax, err);
    %[C(:,i),Id(:,i),~] =Nearest_Subspace_Neighbor(Xi, X(:,i), kmax, size(X,1), i);
end
%% recover real C
% for i=1:1:size(C,2)
%     for j=1:1:size(C,1)
%         if (Id(j,i)~=0)
%             C_spectral_1(i,floor(Id(j,i)))=abs(C(j,i));
%         end
%     end
% end

n=max(s);
%n=3;
N=size(X,2);
% CKSym = BuildAdjacency(thrC(C_spectral_1,rho));
CKSym=abs(C)+abs(C');
grps = SpectralClustering(CKSym,n);
group_real=recover_real_group(N,w,grps);

missrate = Misclassification(group_real,s);

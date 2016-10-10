%% This file read the image, normlize the image
% % collaps the matrix and run SSC_OMP on the data set
% % kmax is the max number of other point we need
% % Xnor is the normalized input data set
% %
%% read the data set
clear;
load('gray_img.mat');
bg=zeros(size(Ig1,1),size(Ig1,2),6);
for i=2:1:7
    bg(:,:,i-1)=eval(['Ig' num2str(i)]);
end
X=reshape(bg,size(Ig1,1)*size(Ig1,2),6)';
load('ground_truth_4.mat');
%% normalized the data set


%% collapse the matrix
load('Matrix_Collapse2.mat');
% ep=7;
% [X,Xs,w]=Matrix_Collapse(X,ep);

%% Sparse subspce clustering
Nx=repmat(sqrt(sum(Xs.*Xs,1)),6,1);
Xsn=Xs./Nx;
kmax=10; err=0.001;
%[C,Id,k]=SparseSubspace(Xsn,kmax,err);
rho = 0.7;
[C,Id,missrate,group_real]=SparseSubspace(Xsn,kmax,err,img,rho,w);

% C_spectral_1=zeros(size(C,2),size(C,2));
% 
% %% recover the big C matrix
% for i=1:1:size(C,2)
%     for j=1:1:size(C,1)
%         C_spectral_1(i,floor(Id(j,i)))=abs(C(j,i));
%     end
% end
% 
% %% spectral clustering
% C_spectral=C_spectral_1'+C_spectral_1;
% groups = SpectralClustering(C_spectral,5);
% 
% %% recover the group_real
% group_real=zeros(size(Xnor,2),1);
% z=zeros(size(w,1),1);
% w_new=[w z];
% for i=1:1:size(groups,1)
%     m=1;
%     while(w_new(i,m)~=0)
%         group_real(w_new(i,m),1)=groups(i,1);
%         m=m+1;
%     end
% end

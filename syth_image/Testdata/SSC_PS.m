%% This file read the image, normlize the image
% % collaps the matrix and run SSC_OMP on the data set
% % kmax is the max number of other point we need
% % Xnor is the normalized input data set
% %
%% read the data set

load('gray_img.mat');
bg=zeros(size(Ig1,1),size(Ig1,2),6);
for i=2:1:7
    bg(:,:,i-1)=eval(['Ig' num2str(i)]);
end
X=reshape(bg,size(Ig1,1)*size(Ig1,2),6)';
load('ground_truth.mat');
%% collapse the matrix
ep=15;
[X,Xs,w]=Matrix_Collapse(X,ep);
Nx=repmat(sqrt(sum(Xs.*Xs,1)),6,1);
Xsn=Xs./Nx;

%% use normal SSC to do the clustering

alpha = 800;
 r = 0; affine = true; outlier = false; rho = 0.7;
 [group_real, missrate1,C1] = SSC(Xsn,r,affine,alpha,outlier,rho,img,w);

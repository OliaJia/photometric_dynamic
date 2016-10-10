function [time,varsub,numOut,num_spa]=ransac_fun(X,thr,Tnum,L,W)
% clearvars -except lt_sel
% Data_Read;

% X=zeros(size(Ig,1)-1,size(Ig,2));
% X=Ig(2:size(Ig,1),:);
% X=im2double(X);
 %clear;
 addpath('E:\Spring2016_Harvard\Research299r\syth_image\Testdata\')
 load('gray_img.mat');
% bg=zeros(size(Ig1,1),size(Ig1,2),7);
% for i=1:1:7
%     bg(:,:,i)=eval(['Ig' num2str(i)]);
% end
% X=reshape(bg,size(Ig1,1)*size(Ig1,2),7)';
tstart=tic;
lt_sel=1:1:size(X,1);
Xnew=X(lt_sel,:);
Out=1:1:size(Xnew,2);
numOut=size(Out,2);
subspace=zeros(100,size(Xnew,2));
s=1;
R=100;
% thr=8e0;
% Tnum=3000;
N=size(Xnew,2);
tick=1;
while (tick<8)
    seed=zeros(3,100);
    inliers=zeros(1,100);
    Lmatrix=Xnew(:,Out);
    for it=1:1:R
        P=randperm(size(Out,2));
        seed(:,it)=P(1:3)';
        [l,d,n]=svd([Lmatrix(:,floor(seed(1,it))) Lmatrix(:,floor(seed(2,it))) Lmatrix(:,floor(seed(3,it)))]);
        L3=l*sqrt(d);
        N3=sqrt(d)*n;
        Nhat=L3\Lmatrix;
        Ihat=L3*Nhat;
        Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
        inliers(1,it)=size(find(Res<thr),2);
    end
    [numLinliers,ind]=max(inliers);
    [l,d,n]=svd([Lmatrix(:,floor(seed(1,ind))) Lmatrix(:,floor(seed(2,ind))) Lmatrix(:,floor(seed(3,ind)))]);
    L3=l*sqrt(d);
    N3=sqrt(d)*n;
    Nhat=L3\Lmatrix;
    Ihat=L3*Nhat;
    Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
    Linliers(1,1:numLinliers)=Out(1,find(Res<thr));
    Loutliers(1,1:numOut-numLinliers)=Out(1,find(Res>=thr));
    clear Out;
    subspace(s,1:numLinliers)=Linliers(1,1:numLinliers);
    s=s+1;
    Out=Loutliers;
    numOut=size(Out,2);
    clear Loutliers Linliers Lmatrix
    thr=thr*3;
    numin(tick,1)=numLinliers;
    tick=tick+1;
end
subspace(s,1:size(Out,2))=Out;
s=s+1;
num_spa=s;
Final_sub=subspace(1:s+1,:);
label=zeros(size(Xnew,2),1);
for i=1:1:num_spa-1
    for j=1:1:N
        if(Final_sub(i,j)==0)
            break;
        end
        label(floor(Final_sub(i,j)),1)=i;
    end
    vis_sub=uint8(reshape(label,L,W)*20);
    imshow(vis_sub);
    pause(0.5);
end
varsub=var(numin);
time=toc(tstart);

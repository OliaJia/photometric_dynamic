clear;
load('gray_img.mat');
bg=zeros(size(Ig1,1),size(Ig1,2),7);
for i=1:1:7
    bg(:,:,i)=eval(['Ig' num2str(i)]);
end
X=reshape(bg,size(Ig1,1)*size(Ig1,2),7)';
L=size(X,1);
N=size(X,2);
illum=sum(X.*X,2);
normX=X./repmat(sqrt(illum),[1,N]);
lightc=zeros(7,1);
lightv=zeros(1,N);
numl=7;
lightset=zeros(30,numl);
subspace=zeros(30,N);
thr=0.6e-14;
% counter of subspace, it is really important, avoid interupting it
s=1;
%% First select the first three lights
% %get the order to select the light
k=1;
while(k<=7)
    if k==1
        [~,lightc(k,1)]=max(illum);
    else
        a=sum(normX*lightv',2);
        neg=illum;
        c=a.*neg;
        c(c==0)=10^20;
        [cs,i]=sort(c);
        lightc(k,1)=i(floor((7-k)/2)+1,1);
    end
    lightv(k,:)=X(floor(lightc(k,1)),:);
    normX(floor(lightc(k,1)),:)=0;
    X(floor(lightc(k,1)),:)=0;
    k=k+1;
end

%% we begin to find out the
k=3;

seed=zeros(3,100);
inliers=zeros(1,100);
for it=1:1:20
    P=randperm(N);
    seed(:,it)=P(1:3)';
    [l,d,n]=svd([lightv(1:3,floor(seed(1,it))) lightv(1:3,floor(seed(2,it))) lightv(1:3,floor(seed(3,it)))]);
    L3=l*sqrt(d);
    N3=sqrt(d)*n;
    Nhat=L3\lightv(1:3,:);
    Ihat=L3*Nhat;
    Res=sqrt(sum((lightv(1:3,:)-Ihat).*(lightv(1:3,:)-Ihat),1));
    inliers(1,it)=size(find(Res<thr),2);
end
%% Record the inliers
[~,ind]=max(inliers);
[l,d,n]=svd([lightv(1:3,floor(seed(1,ind))) lightv(1:3,floor(seed(2,ind))) lightv(1:3,floor(seed(3,ind)))]);
L3=l*sqrt(d);
N3=sqrt(d)*n;
Nhat=L3\lightv(1:3,:);
Ihat=L3*Nhat;
Res=sqrt(sum((lightv(1:3,:)-Ihat).*(lightv(1:3,:)-Ihat),1));

IN=find(Res<0.5e-14);
subspace(s,1:size(IN,2))=IN;
lightset(s,1:3)=lightc(1:3,1)';
s=s+1;
OUT=find(Res>=thr);
%IN3lightv=lightv(1:3,IN);
%OUT3lightv=lightv(1:3,OUT);
numOUT=size(OUT,2);  %%this is an important variable
%% begin to process the inliers
%     masklight=lightv(4:7,IN);
%     masklight(masklight>0)=1;
%     subspace_num=zeros(1,size(masklight,2));
%     for e=0:1:3
%         subspace_num=subspace_num+masklight(e+1,:).*2^e;
%     end
%     subspace=zeros(30,70000);
%     for sub=0:1:15
%         a=find(subspace_num==sub);
%         subspace(sub+1,1:size(a,2))=a;
%     end
%
%% add a light
lselect=zeros(numl-3,1);
for i=1:1:4
    pos=size(find(lightv(i+3,IN)>0),2);
    neg=size(find(lightv(i+3,IN)==0),2);
    lselect(i,1)=abs(pos-neg);
end
[~,indlight]=min(lselect);
lightcnew(k+1,1)=lightc(indlight+3,1);
sub1= find(lightv(indlight+3,IN)>0);
subspace(s,1:size(sub1,2))=sub1;
lightset(s,1:k+1)=lightc([1:3 indlight+3],1);

s=s+1;
sub2= find(lightv(indlight+3,IN)==0);
subspace(s,1:size(sub2,2))=sub2;
lightset(s,1:k)=lightc(1:3,1);
s=s+1;
%% Try this light on outliers
lightry=[1,2,indlight+3; 1,3,indlight+3; 2,3,indlight+3];
numLinliers=zeros(3,1);
Linliers=zeros(3,size(OUT,2));
Loutliers=zeros(3,size(OUT,2));
Lind=zeros(3,1);

for tr=1:1:3
    seed=zeros(3,100);
    inliers=zeros(1,100);
    for it=1:1:20
        Lmatrix=lightv(lightry(tr,:),OUT);
        P=randperm(size(OUT,2));
        seed(:,it)=P(1:3)';
        [l,d,n]=svd([Lmatrix(:,floor(seed(1,it))) Lmatrix(:,floor(seed(2,it))) Lmatrix(:,floor(seed(3,it)))]);
        L3=l*sqrt(d);
        N3=sqrt(d)*n;
        Nhat=L3\Lmatrix;
        Ihat=L3*Nhat;
        Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
        inliers(1,it)=size(find(Res<thr),2);
    end
    [numLinliers(tr,1),ind(tr,1)]=max(inliers);
    [l,d,n]=svd([Lmatrix(:,floor(seed(1,ind(tr,1)))) Lmatrix(:,floor(seed(2,ind(tr,1)))) Lmatrix(:,floor(seed(3,ind(tr,1))))]);
    L3=l*sqrt(d);
    N3=sqrt(d)*n;
    Nhat=L3\Lmatrix;
    Ihat=L3*Nhat;
    Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
    Linliers(tr,1:numLinliers(tr,1))=OUT(1,find(Res<thr));
    Loutliers(tr,1:numOUT-numLinliers(tr,1))=OUT(1,find(Res>=thr));
end
[~,winner]=max(numLinliers);
lightset(s,1:size(lightry,2))=lightc(lightry(winner,:),1);
subspace(s,1:numLinliers(winner,1))=Linliers(winner,1:numLinliers(winner,1));
s=s+1;
clear OUT;
OUT=Loutliers(winner,1:numOUT-numLinliers(winner,1));
numOUT=size(OUT,2);

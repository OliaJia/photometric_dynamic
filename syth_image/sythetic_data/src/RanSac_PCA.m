function [L_est,N_est,IN,OUT]=RanSac_PCA(I,R)
num_sample=4;
Thr=0.3*sum(mad(I'))/size(I,1);
seed=zeros(num_sample,R);
inliers=0;
num=size(I,2);
for it=1:1:R
    P=randperm(num);
    seed(:,it)=P(1:num_sample)';
    [l,d,n]=svd([I(:,floor(seed(1,it)))...
        I(:,floor(seed(2,it))) I(:,floor(seed(3,it)))]);
    L=l*sqrt(d);
    N=sqrt(d)*n;
    Nhat=L\I;
    Ihat=L*Nhat;
    Res=sqrt(sum((I-Ihat).*(I-Ihat),1));
    if(size(find(Res<=Thr),2)>inliers)
        inliers=size(find(Res<=Thr),2);
        IN=find(Res<Thr);
        OUT=find(Res>Thr);
        N_est=Nhat(:,IN);
        L_est=L;
    end
end

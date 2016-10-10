% ransac for regression
% Input:
% I: input matrix (Y in regression problem);
% R: iteration of Ransac;
% Norm: known normal of pixels (X in regression problem);
% Thr: threshold for ransac;
% Output:
% L_est: estimated light direction
% IN: the index of pixels that belongs to the major model


function [L_est,IN,Res]=RanSac_Reg(I,R,Norm)
Thr=0.3*sum(mad(I'))/size(I,1);
seed=zeros(3,R);
inliers=0;
num=size(I,2);
for it=1:1:R
    P=randperm(num);
    seed(:,it)=P(1:3)';
    
    I_sample=[I(:,floor(seed(1,it)))...
        I(:,floor(seed(2,it))) I(:,floor(seed(3,it)))];
    
    N_sample=[Norm(:,floor(seed(1,it)))...
        Norm(:,floor(seed(2,it))) Norm(:,floor(seed(3,it)))];
    
    Lhat=I_sample/N_sample;
    Ihat=Lhat*Norm;
    Res=sqrt(sum((I-Ihat).*(I-Ihat),1));
    if(size(find(Res<=Thr),2)>inliers)
        inliers=size(find(Res<=Thr),2);
        IN=find(Res<Thr);
    end
end
L_est=Lhat;

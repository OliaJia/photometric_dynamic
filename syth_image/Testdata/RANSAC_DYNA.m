%% load data, generate

% Futher question: what if I change the order of adding light? 
% if I use a new set of data
%
clear;
load('gray_img.mat');
bg=zeros(size(Ig1,1),size(Ig1,2),7);
for i=1:1:7
    bg(:,:,i)=eval(['Ig' num2str(i)]);
end
X=reshape(bg,size(Ig1,1)*size(Ig1,2),7)';
M=size(X,1);
N=size(X,2);
illum=sum(X.*X,2);
normX=X./repmat(sqrt(illum),[1,N]);

%% Initialize the parameter 
% % lightc is is a list of light
% % lightv is the corresponding X
% %lightset is to store the lightset that illuminate a subspace
% % subspace is to store the subspace we get in the process
% % subr is to store the range in subspace of the current subspace result
lightc=zeros(M,1);
lightv=zeros(M,N);
lightset=zeros(100,M);
subspace=zeros(100,N);
subr=ones(2,1);
num_spa=0;
thr=1e-14;
TH=100;
R=100;
s=1; % counter in the subspace, the most important counter in this program
%% Part 1 generate a list of light according to some rule
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

%% Run first Ransac on the situation with first three light
% % first set the dynamic parameter lt_sel,lt_left, which is the selected
% light and the left light
lt_sel=(1:3);
lt_left=(4:M);

seed=zeros(3,100);
inliers=zeros(1,100);
for it=1:1:R
    P=randperm(N);
    seed(:,it)=P(1:3)';
    [l,d,n]=svd([lightv(lt_sel,floor(seed(1,it))) lightv(lt_sel,floor(seed(2,it))) lightv(lt_sel,floor(seed(3,it)))]);
    L3=l*sqrt(d);
    N3=sqrt(d)*n;
    Nhat=L3\lightv(lt_sel,:);
    Ihat=L3*Nhat;
    Res=sqrt(sum((lightv(lt_sel,:)-Ihat).*(lightv(lt_sel,:)-Ihat),1));
    inliers(1,it)=size(find(Res<thr),2);
end
%% Record the inliers
[~,ind]=max(inliers);
[l,d,n]=svd([lightv(lt_sel,floor(seed(1,ind))) lightv(lt_sel,floor(seed(2,ind))) lightv(lt_sel,floor(seed(3,ind)))]);
L3=l*sqrt(d);
N3=sqrt(d)*n;
Nhat=L3\lightv(lt_sel,:);
Ihat=L3*Nhat;
Res=sqrt(sum((lightv(lt_sel,:)-Ihat).*(lightv(lt_sel,:)-Ihat),1));

IN=find(Res<thr);
subr(1,1)=s;
%note lightset is updated with subspace always
subspace(s,1:size(IN,2))=IN;
lightset(s,1:size(lt_sel,2))=lightc(lt_sel,1)';
s=s+1;
subr(2,1)=s-1;
% record the number of space
num_spa=subr(2,1)-subr(1,1)+1;
OUT=find(Res>=thr);
numOUT=size(OUT,2);  %%this is an important variable

%% Part 3 begin the interation: in each iteration do: select new light; save point in OUT
while(size(lt_sel,2)<8 && numOUT>TH )
    %% select the next light
    % % lselect is a space helping select light sources
    lselect=zeros(2*num_spa,size(lt_left,2));
    for i=1:1:size(lt_left,2)
        for j=1:1:num_spa
            temp=subspace(subr(1,1)+j-1,:);
            in1=temp(temp>0);
            lselect((j-1)*2+1,i)=size(find(lightv(lt_left(i),in1)>0),2);
            temp=subspace(subr(1,1)+j-1,:);
            in2=temp(temp>0);
            lselect((j-1)*2+2,i)=size(find(lightv(lt_left(i),in2)==0),2);
        end
    end
    clear temp in1 in2;
    [~,win]=min(std(lselect));
    sub_pre=subr;
    subr(1,1)=s;
    for j=1:1:num_spa
        temp=subspace(sub_pre(1,1)+j-1,:);
        in1=temp(temp>0);
        mask1=zeros(1,N); mask1(1,in1)=1;
        mask11=lightv(lt_left(win),:); mask11(mask11>0)=1;
        sub1=find(mask1.*mask11>0); %% have problem
        subspace(s,1:size(sub1,2))=sub1;
        templ=lightset(sub_pre(1,1)+j-1,:);
        inl1=templ(templ>0); inlnew1=[inl1 lightc(lt_left(win),1)];
        lightset(s,1:size(inlnew1,2))=inlnew1;
        s=s+1;
        temp=subspace(sub_pre(1,1)+j-1,:);
        in1=temp(temp>0);
        mask2=zeros(1,N); mask2(1,in1)=1;
        mask22=ones(1,N)-mask11; 
        sub2=find(mask2.*mask22>0); %% have problem
        subspace(s,1:size(sub2,2))=sub2;
        templ=lightset(sub_pre(1,1)+j-1,:);
        inl2=templ(templ>0); inlnew2=inl2;
        lightset(s,1:size(inlnew2,2))=inlnew2;
        s=s+1;
    end
    clear temp templ in1 in2 inl1 inl2 inlnew1 inlnew2;
    subr(2,1)=s-1;
    num_spa=subr(2,1)-subr(1,1)+1;
    ns=size(lt_sel,2); lt_sel(1,ns+1)=lt_left(1,win);
    lt_left(1,win)=0;
    temp=lt_left(lt_left>0);
    clear lt_left; lt_left=temp; clear temp;
    
    %% now use the new light to deal with the issure of OUT
    % % first select 3 possible comination
    com=3;
    lightry=zeros(com,3);
    for i=1:1:com
        p=randperm(size(lt_sel,2)-1);
        lightry(i,:)=[p(1:2) lt_sel(1,size(lt_sel,2))];
    end
    numLinliers=zeros(com,1);
    Linliers=zeros(com,size(OUT,2));
    Loutliers=zeros(com,size(OUT,2));
    Lind=zeros(com,1);
    %for every selected set of light do ransac to see the number of inliers
    for tr=1:1:com
        seed=zeros(3,100);
        inliers=zeros(1,100);
        for it=1:1:R
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
    subr(2,1)=s-1;
    num_spa=subr(2,1)-subr(1,1)+1;
    clear OUT;
    OUT=Loutliers(winner,1:numOUT-numLinliers(winner,1));
    numOUT=size(OUT,2);
end
% record the subspace and display them
subspace(s,1:size(OUT,2))=OUT;
s=s+1;
subr(2,1)=s-1;
num_spa=subr(2,1)-subr(1,1)+1;
Final_sub=subspace(subr(1,1):subr(2,1),:);
label=zeros(N,1);
for i=1:1:num_spa
    for j=1:1:N
        if(Final_sub(i,j)==0)
            break;
        end
        label(floor(Final_sub(i,j)),1)=i;
    end
    vis_sub=uint8(reshape(label,size(Ig1,1),size(Ig1,2))*20);
    imshow(vis_sub);
    pause(0.5);
end





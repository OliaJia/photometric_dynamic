%% load data, generate

% Futher question: what if I change the order of adding light?
% if I use a new set of data
%
clear;
% load('gray_img.mat');
% bg=zeros(size(Ig1,1),size(Ig1,2),7);
% for i=1:1:7
%     bg(:,:,i)=eval(['Ig' num2str(i)]);
% end
% X=reshape(bg,size(Ig1,1)*size(Ig1,2),7)';
clear;
addpath('E:\Spring2016_Harvard\Research299r\syth_image\Testdata\')
load('gray_img.mat');
bg=zeros(size(Ig1,1),size(Ig1,2),7);
for i=1:1:7
    bg(:,:,i)=eval(['Ig' num2str(i)]);
end
X=reshape(bg,size(Ig1,1)*size(Ig1,2),7)';
M=size(X,1);
N=size(X,2);
illum=sum(X.*X,2);
%normX=X./repmat(sqrt(illum),[1,N]);
normX=X;
%% Initialize the parameter
% % lightc is is a list of light
% % lightv is the corresponding X
% %lightset is to store the lightset that illuminate a subspace
% % subspace is to store the subspace we get in the process
% % subr is to store the range in subspace of the current subspace result
%lightc=zeros(M,1);
% actualy this is X
lightv=normX;

lightset=zeros(100,M);
subspace=zeros(100,N);
subr=ones(2,1);
num_spa=0;
thr=1e-19;
TH=500;
R=100;
s=1; % counter in the subspace, the most important counter in this program

%%  Part1: Initial lights selection with matrix column selection
lt_left=1:1:M;
[lt_sel]=LightSel(normX', 3)';
lt_left(lt_sel)=[];
%% Run first Ransac on the situation with first three light
% % first set the dynamic parameter lt_sel,lt_left, which is the selected
% light and the left light
%lt_sel=(1:3);


seed=zeros(3,100);
inliers=zeros(1,100);
for it=1:1:R
    P=randperm(N);
    seed(:,it)=P(1:3)';
    [l,d,n]=svd([lightv(lt_sel,floor(seed(1,it))) lightv(lt_sel,floor(seed(2,it))) lightv(lt_sel,floor(seed(3,it)))]);
    L3=l*sqrt(d);
    if rank(L3)<3
        inliers(1,it)=0;
    else 
        N3=sqrt(d)*n;
        Nhat=L3\lightv(lt_sel,:);
        Ihat=L3*Nhat;
        Res=sqrt(sum((lightv(lt_sel,:)-Ihat).*(lightv(lt_sel,:)-Ihat),1));
        inliers(1,it)=size(find(Res<=thr),2);
    end
end
%% Record the inliers
[~,ind]=max(inliers);
[l,d,n]=svd([lightv(lt_sel,floor(seed(1,ind))) lightv(lt_sel,floor(seed(2,ind))) lightv(lt_sel,floor(seed(3,ind)))]);
L3=l*sqrt(d);
N3=sqrt(d)*n;
Nhat=L3\lightv(lt_sel,:);
Ihat=L3*Nhat;
Res=sqrt(sum((lightv(lt_sel,:)-Ihat).*(lightv(lt_sel,:)-Ihat),1));

IN=find(Res<=thr);
subr(1,1)=s;
%note lightset is updated with subspace always
subspace(s,1:size(IN,2))=IN;
lightset(s,1:size(lt_sel,2))=lt_sel';
s=s+1;
subr(2,1)=s-1;
% record the number of space
num_spa=subr(2,1)-subr(1,1)+1;
OUT=find(Res>thr);
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
            lselect((j-1)*2+1,i)=size(find(lightv(lt_left(i),in1)>0),2);  % cut the existing subspace, which can be improved
            temp=subspace(subr(1,1)+j-1,:);
            in2=temp(temp>0);
            lselect((j-1)*2+2,i)=size(find(lightv(lt_left(i),in2)==0),2);
        end
    end
    clear temp in1 in2;
    score_in=std(lselect);
    % this is the score of light select got from inliers
    %% the fact of outliers to select the light
    seed=zeros(3,100);
    inliers=zeros(1,100);
    lnum=zeros(size(X,1),100);
    for it=1:1:R
        Lmatrix=lightv(:,OUT);
        P=randperm(size(OUT,2));
        seed(:,it)=P(1:3)';
        [l,d,n]=svd([Lmatrix(:,floor(seed(1,it))) Lmatrix(:,floor(seed(2,it))) Lmatrix(:,floor(seed(3,it)))]);
        L3=l*sqrt(d);
        L3_2=abs(L3);
        L3num=sum(L3_2,2);
        lnum(:,it)=(L3num>1e-10);
        N3=sqrt(d)*n;
        Nhat=L3\Lmatrix;
        Ihat=L3*Nhat;
        Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
        inliers(1,it)=size(find(Res<thr),2);
    end
    %%the fact of outliers to select the light
    [index]=LightSel(lightv(:,OUT)', size(lightv,1));
    score_out=zeros(1,size(score_in,2));
    for lrank=1:1:size(lt_left,2)
        score_out(lrank)=find(index==lt_left(lrank));
    end
    score_out=score_out.^2*5;      %score of light selection got from outliers
    [~,win]=min(score_in+score_out);
    
    
    sub_pre=subr;
    subr(1,1)=s;
    % generate new subspaces
    for j=1:1:num_spa
        temp=subspace(sub_pre(1,1)+j-1,:);
        in1=temp(temp>0);
        mask1=zeros(1,N); mask1(1,in1)=1;
        mask11=lightv(lt_left(win),:); mask11(mask11>0)=1;
        sub1=find(mask1.*mask11>0); %% have problem
        subspace(s,1:size(sub1,2))=sub1;
        templ=lightset(sub_pre(1,1)+j-1,:);
        inl1=templ(templ>0); inlnew1=[inl1 lt_left(win)];
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
    % add new light the the lt_sel and delete from lt_left
    ns=size(lt_sel,2); lt_sel(1,ns+1)=lt_left(1,win);
    lt_left(1,win)=0;
    temp=lt_left(lt_left>0);
    clear lt_left; lt_left=temp; clear temp;
    
    %% now use the new light to deal with the issure of OUT
    % % first select 3 possible comination
    %     com=3;
    %     lightry=zeros(com,3);
    %     for i=1:1:com
    %         p=randperm(size(lt_sel,2)-1);
    %         lightry(i,:)=[lt_sel(p(1:2)) lt_sel(1,size(lt_sel,2))];
    %     end
    %     numLinliers=zeros(com,1);
    %     Linliers=zeros(com,size(OUT,2));
    %     Loutliers=zeros(com,size(OUT,2));
    %     Lind=zeros(com,1);
    %     %for every selected set of light do ransac to see the number of inliers
    %     for tr=1:1:com
    %         seed=zeros(3,100);
    %         inliers=zeros(1,100);
    %         for it=1:1:R
    %             Lmatrix=lightv(lightry(tr,:),OUT);
    %             P=randperm(size(OUT,2));
    %             seed(:,it)=P(1:3)';
    %             [l,d,n]=svd([Lmatrix(:,floor(seed(1,it))) Lmatrix(:,floor(seed(2,it))) Lmatrix(:,floor(seed(3,it)))]);
    %             L3=l*sqrt(d);
    %             N3=sqrt(d)*n;
    %             Nhat=L3\Lmatrix;
    %             Ihat=L3*Nhat;
    %             Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
    %             inliers(1,it)=size(find(Res<thr),2);
    %         end
    %         [numLinliers(tr,1),ind(tr,1)]=max(inliers);
    %         [l,d,n]=svd([Lmatrix(:,floor(seed(1,ind(tr,1)))) Lmatrix(:,floor(seed(2,ind(tr,1)))) Lmatrix(:,floor(seed(3,ind(tr,1))))]);
    %         L3=l*sqrt(d);
    %         N3=sqrt(d)*n;
    %         Nhat=L3\Lmatrix;
    %         Ihat=L3*Nhat;
    %         Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
    %         Linliers(tr,1:numLinliers(tr,1))=OUT(1,find(Res<thr));
    %         Loutliers(tr,1:numOUT-numLinliers(tr,1))=OUT(1,find(Res>=thr));
    %     end
    %     [~,winner]=max(numLinliers);
    %     lightset(s,1:size(lightry,2))=lightry(winner,:);
    %% use column selection to choose the three light in the oulier
    
    % find the column selection result from the selected light
    [index2]=ReSelectLight(lightv(lt_sel,OUT)', 2, size(lt_sel,2));
    
    
    %% run ransac to see  the number of inliers in the  selected three lights
    
    saverlight=[lt_sel(1,index2) lt_sel(length(lt_sel))];
    
    seed=zeros(3,100);
    inliers=zeros(1,100);
    for it=1:1:R
        Lmatrix=lightv(saverlight,OUT);
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
    [numLinliers,ind]=max(inliers);
    [l,d,n]=svd([Lmatrix(:,floor(seed(1,ind))) Lmatrix(:,floor(seed(2,ind))) Lmatrix(:,floor(seed(3,ind)))]);
    L3=l*sqrt(d);
    N3=sqrt(d)*n;
    Nhat=L3\Lmatrix;
    Ihat=L3*Nhat;
    Res=sqrt(sum((Lmatrix-Ihat).*(Lmatrix-Ihat),1));
    Linliers(1,1:numLinliers)=OUT(1,find(Res<thr));
    Loutliers(1,1:numOUT-numLinliers)=OUT(1,find(Res>=thr));
    lightset(s,1:size(saverlight,2))=saverlight;
    subspace(s,1:numLinliers)=Linliers(1,1:numLinliers);
    s=s+1;
    subr(2,1)=s-1;
    num_spa=subr(2,1)-subr(1,1)+1;
    clear OUT;
    OUT=Loutliers(1,1:numOUT-numLinliers);
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
    vis_sub=uint8(reshape(label,size(Ig1,1),size(Ig1,2))*10);
    imshow(vis_sub);
    pause(0.5);
end






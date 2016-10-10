%% load data, generate

% Futher question: what if I change the order of adding light?
% if I use a new set of data

%clear;
% Data_Read;
Read_Vis;
X=normc(Ig')';
%X=im2double(X);
M=size(X,1);
N=size(X,2);

% normalize the input image matrix
% normX=X./repmat(sqrt(illum),[1,N]);
normX=X;


%% Initialize the parameter
% % lightc is a list of light
% % lightv is the corresponding X
% % lightset is to store the lightset that illuminate a subspace
% % subspace is to store the subspace we get in the process
% % subr is to store the range in subspace of the current subspace result
%   lightc=zeros(M,1);
% actualy this is X


lightv=normX;
N_estimate=zeros(3,N);
L_estimate=zeros(M,3);
lightset=zeros(100,M);
subspace=zeros(100,N);
subr=ones(2,1);
num_spa=0;
thr=1.3e-17;
thr_light=0.025;
TH=500;
R=100;
s=1; % counter in the subspace, the most important counter in this program
alpha=0.002;
beta=0.5;
%%  Part1: Initial lights selection with matrix column selection
lt_left=1:1:M;
[lt_sel]=LightSel(normX', 4)';
lt_left(lt_sel)=[];
%% Run first Ransac on the situation with first three light
% % first set the dynamic parameter lt_sel,lt_left, which is the selected
% light and the left light

I=lightv(lt_sel,:);
[L_est,N_est,IN,OUT]=RanSac_PCA(lightv(lt_sel,:),R);

light_norm=sqrt(sum(L_est.*L_est,2));
lt_vis=lt_sel(light_norm>thr_light);
N_estimate(:,IN)=N_est;
L_estimate(lt_vis,:)=L_est;


subr(1,1)=s;

%note lightset is updated with subspace always
subspace(s,1:size(IN,2))=IN;



lightset(s,1:size(lt_sel,2))=lt_vis';
s=s+1;
subr(2,1)=s-1;

% record the number of space
num_spa=subr(2,1)-subr(1,1)+1;

numOUT=size(OUT,2);  %%this is an important variable

%% Part 3 begin the interation: in each iteration do: select new light; save point in OUT
while(size(lt_sel,2)<8 && numOUT>TH)         
    lselect=zeros(2*num_spa,size(lt_left,2));
    Res_records=zeros(num_spa,size(IN,2),size(lt_left,2));
    thr_in=0.002;
    vis_flag=zeros(1,size(lt_sel,2));
    for i=1:1:size(lt_left,2)
        for j=1:1:num_spa
            temp_l=lightset(subr(1,1)+j-1,:);
            elight=temp_l(temp_l>0);
            tempp=subspace(subr(1,1)+j-1,:);
            epixel=tempp(tempp>0);
            temp_lightv=lightv([lt_sel,lt_left(i)],epixel);
            temp_N=N_estimate(:,epixel);
            % start Ransac to cut the subspace
            [L_est,inliers_temp,Res_records(j,:,i)]=RanSac_Reg(temp_lightv,R,temp_N);
            vis_flag(1,i)=(norm(L_est(size(L_est,1),:),2)>thr_light);
            inliers_temp_num=size(inliers_temp,2);
            lselect((j-1)*2+1,i)=inliers_temp_num;
            lselect((j-1)*2+2,i)=size(epixel,2)-lselect((j-1)*2+1,i);
        end
    end
    
    %% cut each existing subspaces into two
    Score_in=std(lselect);
    [~,rank]=sort(Score_in,'ascend');
    
    [Score_out,best_in]=light_increment(lightv(:,OUT)',lt_sel,lt_left,3);
    
    Score_select=alpha*Score_in+beta*Score_out';
    
    [~,win_index]=min(Score_select);
    
    win=lt_left(1,win_index);
    
    lt_sel(1,size(lt_sel,2)+1)=lt_left(1,win);
    
    lt_left(win_index)=[];
    % separate the subspaces in last step
    
    

    %% the fact of outliers to select the light    
    sub_pre=subr;
    subr(1,1)=s;
    % generate new subspaces
    for j=1:1:num_spa
        Res_win=Res_records(j,:,win);
        sub1=find(Res_win<thr_in);
        subspace(s,1:size(sub1,2))=sub1;
        lightset(s,:)=lightset(sub_pre(1,1)-1+j,:);
        if vis_flag==1
            lightset(s,size(lightset(s,lightset(s,:)>0),2)+1)=win;
        end
        s=s+1;
        sub2=find(Res_win>=thr_in);
        subspace(s,1:size(sub2,2))=sub2;
        lightset(s,:)=lightset(sub_pre(1,1)-1+j,:);
        if vis_flag==0
            lightset(s,size(lightset(s,lightset(s,:)>0),2)+1)=win;
        end
        s=s+1;
    end
    subr(2,1)=s-1;
    num_spa=subr(2,1)-subr(1,1)+1;
    
    %% use column selection to choose the three light in the oulier
    
    % find the column selection result from the selected light
    [L_est,N_est,IN,OUT]=RanSac_PCA(lightv(lt_sel,OUT),R);
    
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
    vis_sub=uint8(reshape(label,L,W)*20);
    imshow(vis_sub);
    pause(0.5);
end






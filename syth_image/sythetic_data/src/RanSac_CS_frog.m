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

subr=ones(2,1);
num_spa=0;
thr_light=0.025;
TH=500;
R=100;
Thr_Num=5000;
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
OUT=1:1:size(lightv,2);
N_normal=zeros(3,size(lightv,2));

[N_normal, OUT, subspace_new, lightset_new,num_spa]=subspace_discover(OUT,lightv,lt_sel,Thr_Num,N_normal, num_spa);


subspace=subspace_new;
lightset=lightset_new;
while(size(lt_sel,2)<8 && size(OUT,2)>TH)         
    lselect=zeros(2*num_spa,size(lt_left,2));
    Res_records=cell(num_spa,size(lt_left,2));
    thr_in=0.002;
    vis_flag=zeros(num_spa,size(lt_left,2));
    for i=1:1:size(lt_left,2)
        for j=1:1:num_spa
            temp_lightv=lightv([lt_sel,lt_left(i)],subspace{j});
            temp_N=N_normal(:,subspace{j});
            % start Ransac to cut the subspace
            [L_est,inliers_temp,Res_records{i,j}]=RanSac_Reg(temp_lightv,R,temp_N);
            if(norm(L_est(size(L_est,1),:),2)>thr_light)
                vis_flag(j,i)=1;
            else
                vis_flag(j,i)=0;
            end
            inliers_temp_num=size(inliers_temp,2);
            lselect((j-1)*2+1,i)=inliers_temp_num;
            lselect((j-1)*2+2,i)=size(subspace{j},2)-lselect((j-1)*2+1,i);
        end
    end
    Score_in=std(lselect);
    [~,rank]=sort(Score_in,'ascend');
    
    [Score_out,best_in]=light_increment(lightv(:,OUT)',lt_sel,lt_left,3);
    
    Score_select=alpha*Score_in+beta*Score_out';
    
    [~,win_index]=min(Score_select);
    
    flag_light=vis_flag(:,win_index);
    win=lt_left(1,win_index);
    
    lt_sel(1,size(lt_sel,2)+1)=win;
    
    lt_left(win_index)=[];
    
    
    %% cut subspace
  
    lightset_temp=lightset;
    for i=1:1:num_spa
        sub1=find(Res_records{win_index,i}<thr_in);
        sub2=find(Res_records{win_index,i}>thr_in);
        lightset{2*i-1}=[lightset_temp{i} win];
        lightset{2*i}=lightset_temp{i};
        if(flag_light(i,1)==1)
            subspace{2*i-1}=sub1;
            subspace{2*i}=sub2;
        else
            subspace{2*i-1}=sub2;
            subspace{2*i}=sub1;
        end
    end
    num_spa=num_spa*2;
    
    Thr_Num=Thr_Num*0.1;
    [N_normal, OUT, subspace_new, lightset_new,num_spa_new]=subspace_discover(OUT,lightv,lt_sel,Thr_Num,N_normal, 0);
    for j=1:size(subspace_new)
        subspace{num_spa+j}=subspace_new{j};
        lightset{num_spa+j}=lightset_new{j};
    end
end



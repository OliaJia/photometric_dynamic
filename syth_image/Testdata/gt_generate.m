%*** Yaoguang Jia ***
%*** generate ground truth of visibility subspace
%% read data
load('binary_img.mat');
bg=zeros(size(Ib1,1),size(Ib1,2),7);
for i=1:1:7
    bg(:,:,i)=eval(['Ib' num2str(i)]);
end

code_point=reshape(bg,size(Ib1,1)*size(Ib1,2),7);


Y=1/2;
%scatter(1:size(Ib1,1)*size(Ib2,2),code_point);
%in_var=zeros(5,1);
%elbow method to find out the number of clusters
% for k=2:1:5
%     [idx,C,sumd,D]=kmeans(code_point,k);
%     D=D.^2;
%     for k_r=1:1:k
%         in_var(k,1)=in_var(k,1)+sum(D(find(idx==k_r),k_r))/size(find(idx==k_r),1);
%     end
% end
%% cluster the ground truth just collapse the same element (5 subspace)
% % just collapse the matrix with greedy search
ind=1;
for i=1:1:size(code_point,1)
        if  sum(code_point(i,:),2) ~= 0
            w(ind,1) =i;
            p = 2;
            for j = i+1 : 1 : size(code_point,1)
                if sum(code_point(j,:),2) == 0
                    continue;
                else
                    t=sum(abs(code_point(i,:)-code_point(j,:)),2);
                    if  t<=1
                        code_point(j,:)=0;
                        w(ind,p)=j;
                        p=p+1;
                    end
                end
            end
            ind=ind + 1;
        end
end
img=zeros(size(code_point,1),1);

%% display the ground truth of 

for i=1:1:size(w,1)
    for j=1:1:size(w,2)
        if w(i,j)==0
            continue;
        else
            img(floor(w(i,j)),1)=i;
        end
    end
end
vis_sub=uint8(reshape(img,size(Ib1,1),size(Ib1,2))*20);
imshow(vis_sub);
  
% [idx,C,sumd,D]=kmeans(code_point,5);
% vis_sub=uint8(reshape(idx,size(Ib1,1),size(Ib1,2))*40);
% imshow(vis_sub);


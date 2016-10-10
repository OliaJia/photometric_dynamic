%% read data
load('binary_img.mat');
bg=zeros(size(Ib1,1),size(Ib1,2),7);
for i=1:1:7
    bg(:,:,i)=eval(['Ib' num2str(i)]);
end

code_point=reshape(bg,size(Ib1,1)*size(Ib1,2),7);
code_point_2=zeros(size(code_point,1),1);
for i=1:1:7
    code_point_2=code_point_2+code_point(:,i)*2^(i-1);
end
[img,C,sumd,D]=kmeans(code_point_2,4);
vis_sub=uint8(reshape(img,size(Ib1,1),size(Ib1,2))*40);
imshow(vis_sub);

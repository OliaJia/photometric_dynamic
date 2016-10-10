img=zeros(size(X,2),1);
for i=1:1:size(w,1)
    for j=1:1:size(w,2)
        if w(i,j)==0
            continue;
        else
            img(floor(w(i,j)),1)=i;
        end
    end
end
 vis_sub=uint8(reshape(img,size(Ig1,1),size(Ig1,2))*0.8);
imshow(vis_sub);
  
    

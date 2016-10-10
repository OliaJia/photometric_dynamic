

clear;


addpath 'E:\Spring2016_Harvard\Research299r\syth_image\sythetic_data\image2\';
image='E:\Spring2016_Harvard\Research299r\syth_image\sythetic_data\image2\*.png';
imagepath='E:\Spring2016_Harvard\Research299r\syth_image\sythetic_data\image2\';
file=dir(image);
for i=1:1:length(file)
    Ib=imread([imagepath,file(i).name]);
    M=imresize(rgb2gray(Ib),0.6);
    imshow(M);
    Ig(i,:)=reshape(M,1,size(M,1)*size(M,2));
end
Ref1=imread('IMG_000.png');
Ref2=imresize(rgb2gray(Ref1),0.6);
Ref3=reshape(Ref2,1,size(M,1)*size(M,2));

Ig=Ig-repmat(Ref3,length(file),1);
L=size(M,1); W=size(M,2);

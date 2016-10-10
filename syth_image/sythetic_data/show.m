



clear;

addpath 'E:\Spring2016_Harvard\Research299r\syth_image\sythetic_data\image2\';
image='E:\Spring2016_Harvard\Research299r\syth_image\sythetic_data\image2\*.png';
imagepath='E:\Spring2016_Harvard\Research299r\syth_image\sythetic_data\image2\';
file=dir(image);
for i=1:1:length(file)
    Ib=imread([imagepath,file(i).name]);
    %M=imresize(rgb2gray(Ib),0.6)
    %imshow(M);
    Ref1=imread('IMG_000.png');
    
    k=Ib-Ref1;
    figure;
    imshow(k);
end


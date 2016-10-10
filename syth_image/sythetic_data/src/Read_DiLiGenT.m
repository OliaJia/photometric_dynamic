


clear;
bitDepth=16;
gamma=1;

addpath '/Users/Yaoguang/Dropbox/photometric_dynamic/DiLiGenT/pmsData/catPNG/';
image='/Users/Yaoguang/Dropbox/photometric_dynamic/DiLiGenT/pmsData/catPNG/*.png';
imagepath='/Users/Yaoguang/Dropbox/photometric_dynamic/DiLiGenT/pmsData/catPNG/';
file=dir(image);
load 'light_intensities.txt';
mask=imread('/Users/Yaoguang/Dropbox/photometric_dynamic/DiLiGenT/pmsData/catPNG/true_value/mask.png');
m=find(mask>0);

for i=1:1:length(file)
    Ib=imread([imagepath,file(i).name]);
    E = (double(Ib)./ (2^bitDepth-1)).^gamma; % for float input, set bitDepth = 1;
    %E = imresize(E, resize, 'nearest');
	
	[H W C] = size(E);
    % Normalize the image with light intensities
	E = reshape(reshape(E,[H*W C])*diag(1./light_intensities(i,:)),[H W C]);
	E = max(0,E);
    M=rgb2gray(E);
    M=M(m);
    %M=imresize(rgb2gray(Ib),0.6);
    %imshow(M);
    Ig(i,:)=reshape(M,1,size(M,1)*size(M,2));
end
% Ref1=imread('IMG_000.png');
% Ref2=imresize(rgb2gray(Ref1),0.6);
% Ref3=reshape(Ref2,1,size(M,1)*size(M,2));

L=size(E,1); W=size(E,2);

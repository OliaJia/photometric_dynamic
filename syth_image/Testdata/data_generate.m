I0=imread('IMG_000.png');
I1=imread('IMG_001.png');
I2=imread('IMG_002.png');
I3=imread('IMG_003.png');
I4=imread('IMG_004.png');
I5=imread('IMG_005.png');
I6=imread('IMG_006.png');
I7=imread('IMG_007.png');
I1=I1-I0;I2=I2-I0;I3=I3-I0;
I4=I4-I0;I5=I5-I0;I6=I6-I0;
I7=I7-I0;
imshow(I1);
Ig1=rgb2gray(I1);Ig2=rgb2gray(I2);Ig3=rgb2gray(I3);
Ig4=rgb2gray(I4);Ig5=rgb2gray(I5);Ig6=rgb2gray(I6);
Ig7=rgb2gray(I7);
Ig1(find(Ig1>0))=1;Ig2(find(Ig2>0))=1;Ig3(find(Ig3>0))=1;
Ig4(find(Ig4>0))=1;Ig5(find(Ig5>0))=1;Ig6(find(Ig6>0))=1;
Ig7(find(Ig7>0))=1;
Ib1=Ig1;Ib2=Ig2;Ib3=Ig3;Ib4=Ig4;Ib5=Ig5;Ib6=Ig6;Ib7=Ig7;

% This program do the column selection
% with the greedy search algorithm 

%% Generate test data
light=7;
pix=2000;

X=rand(pix,light);

% addpath('E:\Spring2016_Harvard\Research299r\syth_image\dual\dual\')
% 
% lambda=0.4;
% [A,E,~] = dual_rpca(X, lambda );

[U,S,V]=svd(X,0);

UV=U(:,1:3)*S(1:3,1:3)*V(1:3,:);


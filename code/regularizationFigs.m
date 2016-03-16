% demonstrate regularization assumptions

clear all; close all; clc;

cd('..')
addpath(genpath(pwd))
cd('wavelet_FISTA')

% load image
X = double(rgb2gray(imread('peppers.png')));
X = X/255;
[m,n] = size(X);

%% sparsity in wavelet domain

W = opHaar2(m,n,1);
Y = reshape(W*X(:),m,n);
figure(), imshow(Y,[]), colorbar('Location','SouthOutside')

W = opHaar2(m,n,5);
Y = reshape(W*X(:),m,n);
figure(), imshow(Y,[]), colorbar('Location','SouthOutside')

%% low total variation

D1 = abs(diff(X,1,1));
figure(), imshow(D1,[]), colorbar('Location','SouthOutside')

D2 = abs(diff(X,1,2));
figure(), imshow(D2,[]), colorbar('Location','SouthOutside')
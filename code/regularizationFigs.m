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

W = opHaar2(m,n,2);
Y = reshape(W*X(:),m,n);
figure(), imshow(Y,[]), colorbar

%% low total variation

D1 = diff(X,1,1);
D2 = diff(X,1,2);
D = abs(D1(:,1:511) + D2(1:383,:));
figure(), imshow(D,[]), colorbar
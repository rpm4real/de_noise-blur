% demonstrate how choice of lambda affects results

clear all; close all; clc;

cd('..')
addpath(genpath(pwd))
cd('wavelet_FISTA')

% load image
X = double(rgb2gray(imread('peppers.png')));
X = X/255;
[m,n] = size(X);

% blur and add noise
Bpars.nType = 'both';
Bobs = blurt(X,Bpars);

% transforms
W1 = opHaar2(m,n);
W2 = opWavelet2(m,n,'Daubechies');

% fix images
Dpars.fig = 0; Dpars.dispfunc = 0; 
lambdaVals = [1e-4 1e-3 1e-2];
error1 = zeros(1,length(lambdaVals));
error2 = zeros(1,length(lambdaVals));
for i = 1:length(lambdaVals)
    X_out1=deblur_wavelet_2norm(Bobs,W1,lambdaVals(i),1e-4,Dpars);
    X_out2=deblur_wavelet_2norm(Bobs,W2,lambdaVals(i),1e-4,Dpars);
    subplot_tight(2,3,i), imshow(X_out1);
    subplot_tight(2,3,3+i), imshow(X_out2);
    error1(i) = norm(X-X_out1,'fro');
    error2(i) = norm(X-X_out2,'fro');
end
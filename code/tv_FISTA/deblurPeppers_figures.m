% deblur/denoise peppers image using L1 wavelet regularization and FISTA

clear all; close all; clc;

cd('..')
addpath(genpath(pwd))
cd('tv_FISTA')

% load image
X = double(rgb2gray(imread('peppers.png')));
X = X/255;
[m,n] = size(X);

% choose noise: Gaussian = 1, Student's T = 2
noise = 2;

% set wavelet and blur/noise parameters
% WHAT PARAMETERS WERE USED FOR THE FIGURES?
if noise == 1
   Bpars.nType = 'gaussian';
   Bpars.nLevelG = 1e-2;
   lambda1 = 0.003; lambda2 = 0.002; lambda3 = 0.2; gamma1 = 0.02; gamma2 = 50;
else
   Bpars.nType = 'student';
   Bpars.nLevelS = 1e-2;
   lambda1 = 0.05; lambda2 = 0.0025; lambda3 = 0.2; gamma1 = 0.02; gamma2 = 50;
end

% blur and add noise
[P,center] = psfGauss([9,9],1);
Bobs = blurt(X,Bpars);

% deblur/denoise image
Dpars.fig = 0; Dpars.MAXITER = 30; Dpars.denoiseiter = 20;
tic
[X_out1,~] = deblur_tv_fista(Bobs,P,center,lambda1,0,1,Dpars);
toc
norm(X-X_out1,'fro')
tic
[X_out2,~] = HuberTV_Prox(Bobs,gamma1,P,center,lambda2,0,1,Dpars);
toc
norm(X-X_out2,'fro')
tic
[X_out3,~] = LogCoshTVProx(Bobs,gamma2,P,center,lambda3,0,1,Dpars);
toc
norm(X-X_out3,'fro')

% view results (frobenius v. huber)
figure()
subplot_tight(2,2,1), imshow(X)
subplot_tight(2,2,2), imshow(Bobs)
subplot_tight(2,2,3), imshow(X_out1)
subplot_tight(2,2,4), imshow(X_out2)

% view results (logcosh v. huber)
figure()
subplot_tight(2,2,1), imshow(X)
subplot_tight(2,2,2), imshow(Bobs)
subplot_tight(2,2,3), imshow(X_out3)
subplot_tight(2,2,4), imshow(X_out2)
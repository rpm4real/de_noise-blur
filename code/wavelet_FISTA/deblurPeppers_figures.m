% deblur/denoise peppers image using L1 wavelet regularization and FISTA

clear all; close all; clc;

cd('..')
addpath(genpath(pwd))
cd('wavelet_FISTA')

% load image
X = double(rgb2gray(imread('peppers.png')));
X = X/255;
[m,n] = size(X);

% set convergence tolerance
tol = 1e-3;

% choose wavelet and noise
% wav: Haar = 1, Daubechies = 2
% noise: Gaussian = 1, Student's T = 2
wav = 1;
noise = 1;

% set wavelet and blur/noise parameters
if wav == 1
    W = opHaar2(m,n);
    if noise == 1
        Bpars.nType = 'gaussian'; 
        Bpars.nLevelG = 1e-2;
        lambda1 = 1e-3; lambda2 = 1e-3; gamma = 1e-2;
    else
        Bpars.nType = 'student';
        Bpars.nLevelS = 1e-2;
        lambda1 = 1e-1; lambda2 = 1e-2; gamma = 1e-2;
    end
else
    W = opWavelet2(m,n,'Daubechies');
    if noise == 1
        Bpars.nType = 'gaussian'; 
        Bpars.nLevelG = 1e-2;
        lambda1 = 1e-2; lambda2 = 1e-3; gamma = 1e-2;
    else
        Bpars.nType = 'student';
        Bpars.nLevelS = 1e-2;
        lambda1 = 1e-1; lambda2 = 1e-2; gamma = 1e-2;
    end
end

% blur and add noise
Bobs = blurt(X,Bpars);

% deblur/denoise image
Dpars.fig = 0; Dpars.dispfunc = 0; Dpars.MAXITER = 200;
[X_out1] = deblur_wavelet_2norm(Bobs,W,lambda1,tol,Dpars);
[X_out2] = deblur_wavelet_huber(Bobs,W,lambda2,gamma,tol,Dpars);

% view results
figure()
subplot_tight(2,2,1), imshow(X)
subplot_tight(2,2,2), imshow(Bobs)
subplot_tight(2,2,3), imshow(X_out1)
subplot_tight(2,2,4), imshow(X_out2)

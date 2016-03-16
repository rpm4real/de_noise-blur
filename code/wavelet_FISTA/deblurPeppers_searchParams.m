% deblur/denoise peppers image using L1 wavelet regularization and FISTA
% with a search for the best parameters

clear all; close all; clc;

cd('..')
addpath(genpath(pwd))
cd('wavelet_FISTA')

% load image
X = double(rgb2gray(imread('peppers.png')));
X = X/255;
[m,n] = size(X);

% blur and add noise
%Bpars.nType = 'gaussian'; Bpars.nLevelG = 1e-2;
Bpars.nType = 'student'; Bpars.nLevelS = 1e-2;
Bobs = blurt(X,Bpars);

% wavelet transform
%W = opHaar2(m,n);
W = opWavelet2(m,n,'Daubechies');

% tolerance
tol = 1e-3;

% deblurring parameters
Dpars.fig = 0; Dpars.dispfunc = 0; Dpars.MAXITER = 200;

%% find best parameters for 2-norm 

lambdaVals = 10.^(-6:-1);
bestLam = 0; bestError = 100; bestIter = 200;
errors = zeros(1,length(lambdaVals));
for i = 1:length(lambdaVals)
    [X_out,fun_all]=deblur_wavelet_2norm(Bobs,W,lambdaVals(i),tol,Dpars);
    temp = norm(X-X_out,'fro');
    errors(i) = temp;
    if temp < bestError
        bestLam1 = lambdaVals(i);
        bestError = temp;
        bestIter = length(fun_all);
    end
end

%% find best parameters for huber norm

lambdaVals = 10.^(-6:-1);
gammaVals = 10.^(-2:0);
bestLam = 0; bestError = 100; bestIter = 200;
errors = zeros(length(lambdaVals),length(gammaVals));
for i = 1:length(lambdaVals)
    for j = 1:length(gammaVals)
        [X_out,fun_all]=deblur_wavelet_huber(Bobs,W,lambdaVals(i),gammaVals(j),tol,Dpars);
        temp = norm(X-X_out,'fro');
        errors(i,j) = temp;
        if temp < bestError
            bestLam2 = lambdaVals(i);
            bestGam = gammaVals(j);
            bestError = temp;
            bestIter = length(fun_all);
        end
    end
end

%% plot results

[X_out1]=deblur_wavelet_2norm(Bobs,W,bestLam1,tol,Dpars);
[X_out2]=deblur_wavelet_huber(Bobs,W,bestLam2,bestGam,tol,Dpars);

figure()
subplot_tight(2,2,1), imshow(X)
subplot_tight(2,2,2), imshow(Bobs)
subplot_tight(2,2,3), imshow(X_out1)
subplot_tight(2,2,4), imshow(X_out2)


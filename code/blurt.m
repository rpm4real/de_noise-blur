function Bobs = blurt(X,pars)
%--------------------------------------------------------------------------
% This function blurs an image with a Gaussian PSF and reflexive boundary 
% conditions, then adds noise of the specified type and intensity. Input
% can be either a color or black and white, but output is black and white.
%
% INPUT
% 
% X                 True image, either RGB or BW
% pars              Parameter structure
% pars.nType        Noise type: gaussian, student, or both (default = gaussian)
% pars.nLevelG      Gaussian noise level (default = 1e-3)
% pars.nLevelS      Student's T noise level (default = 1e-4)
% pars.bSize        Blur PSF dimension (default = 9)
% pars.bLevel       Blur PSF standard deviation (default = 2)
%
% OUTPUT
% 
% Bobs              Blurred and noisy image (BW)
%
%--------------------------------------------------------------------------

% assign parameters according to pars and/or default values
flag=exist('pars');
if (flag&&isfield(pars,'nType'))
    nType = pars.nType;
else
    nType = 'gaussian';
end 
if (flag&&isfield(pars,'nLevelG'))
    nLevelG = pars.nLevelG;
else
    nLevelG =1e-3;
end
if (flag&&isfield(pars,'nLevelS'))
    nLevelS = pars.nLevelS;
else
    nLevelS = 1e-4;
end
if (flag&&isfield(pars,'bSize'))
    bSize = pars.bSize;
else
    bSize = 9;
end
if (flag&&isfield(pars,'bLevel'))
    bLevel = pars.bLevel;
else
    bLevel = 2;
end

% convert to BW if necessary
if size(X,3) == 3
    X = rgb2gray(X);
end

% scale image if necessary
if max(X(:)) == 255
    X = double(X)/255;
end
[m,n] = size(X);

% blur image
[P,center] = psfGauss([bSize,bSize],bLevel);
A = opConvolve2(m,n,P,center,'reflexive');
b = A*X(:); B = reshape(b,m,n);

% add noise
if ~(strcmp(nType,'gaussian')||strcmp(nType,'student')||strcmp(nType,'both'))
    error('nType should be ''gaussian'', ''student'', or ''both''.') 
end
if (strcmp(nType,'gaussian')||strcmp(nType,'both'))
    B = B + nLevelG*randn(m,n);
end
if (strcmp(nType,'student')||strcmp(nType,'both'))
    B = B + nLevelS*trnd(1,m,n);
    B = B.*(B > 0);
    B = B - (B > 1).*B + (B > 1);
end
Bobs = B;

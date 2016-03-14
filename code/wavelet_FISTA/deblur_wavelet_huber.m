function [X_out,fun_all]=deblur_wavelet_huber(Bobs,W,lambda,gamma,tol,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements FISTA for solving the linear inverse problem with 
% an orthogonal l1 wavelet regularizer and reflexive boundary conditions,
% where the fidelity term uses the Huber norm:
% f(A(X)-Bobs) = min_s 1/2||A(X)-Bobs-s||^2 + gamma||s||_1
% 
% INPUT
%
% Bobs..........................The observed image which is blurred and noisy
% W.............................A spot operator
%                               For an image x = X(:), W*x is an orthogonal
%                               transform of the image X.
%                               W' denotes the inverse transform. 
% lambda........................1-norm regularization parameter
% gamma.........................Huber regularization parameter
% tol...........................Convergence tolerance
% pars..........................Parameters structure
% pars.P........................PSF of the blurring operator 
%                               (default 9x9 gaussian kernel with stdev 2)
% pars.center...................A vector of length 2 containing the center
%                               of the PSF (default [5,5])
% pars.MAXITER..................Maximum number of iterations (Default=100)
% pars.fig......................1 if the image is shown at each iteration, 
%                               0 otherwise (Default=1)
% pars.dispfunc.................1 if you want the function values printed
%                               at each iteration, 0 otherwise (default = 1)
%
% OUTPUT
% 
% X_out ........................Solution of the problem
%                               min{||A(X)-Bobs||^2+lambda \|Wx\|_1
% fun_all ......................Array containing all function values
%                               obtained in the FISTA method
%
%--------------------------------------------------------------------------

% Assigning parameters according to pars and/or default values
flag=exist('pars');
if (flag&&isfield(pars,'P'))
    P=pars.P;
else
    P=psfGauss([9,9],2);
end
if (flag&&isfield(pars,'center'))
    center=pars.center;
else
    center=[5,5];
end
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if (flag&&isfield(pars,'fig'))
    fig=pars.fig;
else
    fig = 1;
end
if (flag&&isfield(pars,'dispfunc'))
    dispfunc=pars.dispfunc;
else
    dispfunc=1;
end

fun_all=[];
[m,n]=size(Bobs);
A=opConvolve2(m,n,P,center,'reflexive');
Pbig=padPSF(P,[m,n]);

% reflexive boundary conditions
trans=@(X)dct2(X);
itrans=@(X)idct2(X);

% compute the eigenvalues of the blurring matrix         
e1=zeros(m,n);
e1(1,1)=1;
Sbig=dct2(dctshift(Pbig,center))./dct2(e1);

% the Lipschitz constant of the gradient of 
% min_s 1/2||A(X)-Bobs - s||^2 + gamma||s||_1
L=max(max(abs(Sbig).^2));

% soft-thresholding
St=@(s)max(abs(s)-gamma,0).*sign(s); % soft-thresholding

% initialization
X_iter=Bobs;
Y=X_iter;
t_new=1;
if dispfunc
    fprintf('************\n');
    fprintf('**FISTA**\n');
    fprintf('************\n');
    fprintf('#iter  fun-val         relative-dif\n==============================\n');
end
for i=1:MAXITER
    % Store the old value of the iterate and the t-constant
    X_old=X_iter;
    t_old=t_new;
   
    % Gradient step
    ry=A*Y(:)-Bobs(:);
    grad=A'*(ry-St(ry));
    Y=Y-1/L*reshape(grad,m,n);
    
    % Wavelet transform
    Wy=W*Y(:);
    % Soft thresholding 
    D=abs(Wy)-lambda/L;
    Wy=sign(Wy).*((D>0).*D);
    % The new iterate inverse wavelet transform of WY
    x_iter=W'*Wy;
    X_iter=reshape(x_iter,m,n);
    
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+(t_old-1)/t_new*(X_iter-X_old);
    
    % Compute the l1 norm of the wavelet transform and the function value and store it in
    % the function values vector fun_all if exists.
    t=sum(sum(abs(W*X_iter(:))));
    rx=A*X_iter(:)-Bobs(:);
    fun_val=1/2*norm(rx-St(rx))+gamma*norm(St(rx),1)+lambda*t;
    fun_all=[fun_all;fun_val];
    % printing the information of the current iteration
    if dispfunc
        fprintf('%3d    %15.5f                %15.5f \n',i,fun_val,norm(X_iter-X_old,'fro')/norm(X_old,'fro'));
    end
        
    if (fig)
        figure(314)
        imshow(X_iter,[])
    end
    
    % check convergence
    if i>1
        if abs(fun_all(i)-fun_all(i-1))/fun_all(i-1)<tol
            i
            break
        end
    end
end

X_out=X_iter;
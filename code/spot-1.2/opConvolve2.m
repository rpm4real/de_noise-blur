classdef opConvolve2 < opSpot
%OPCONVOLVE   Two dimensional convolution operator.

%   opConvolve2(M,N,KERNEL,OFFSET,MODE) creates an operator for
%   two-dimensional convolution. The OFFSET parameter determines the
%   center of the KERNEL and has a default value of [1,1]. There are 
%   four types of MODE:
% 
%   MODE = 'regular'   - convolve input with kernel;
%          'truncated' - convolve input with kernel, but keep only
%                        those MxN entries in the result that
%                        overlap with the input;
%          'periodic'  - do periodic convolution of the input
%          'reflexive' - do reflexive convolution of the input
%
%   The output of the convolution operator, like all other
%   operators, is in vector form.
%
%   Modified (by Kelsey) from opConvolve()
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.
%
%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       funHandle     % Multiplication function
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opConvolve2(m,n,kernel,offset,mode)

          if nargin < 3
             error('opConvolve requires at least three parameters.');
          end

          if nargin < 4, offset = [];  end;
          if nargin < 5, mode = 'regular'; end;
          
          switch lower(mode)
             case {'periodic'}
                periodic = true; reflexive = false; truncated = false;
    
             case {'reflexive'}
                periodic = false; reflexive = true; truncated = false;
             
             case {'truncated'}
                periodic = false; reflexive = false; truncated = true;
                
             case {'regular'}
                periodic = false; reflexive = false; truncated = false;
    
             otherwise
                error('Mode parameter must be one of ''regular'', ''truncated'', ''periodic'', or ''reflexive''.');
          end

          % Get basic information
          if isempty(offset), offset = [1,1]; end;
          if length(offset) < 2
             error('Offset parameter needs to contain two entries.');
          end
          offset = offset(1:2);
          k      = [size(kernel,1),size(kernel,2)];
          cflag  = ~isreal(kernel);
 
          if periodic
             % ========= Periodic convolution =========

             % Zero pad if needed
             bigKernel = zeros(m, n);
             bigKernel(1:size(kernel,1),1:size(kernel,2)) = kernel;
 
             % Precompute kernel in frequency domain
             fKernel = fft2(circshift(bigKernel,1-offset));
 
             % Create function handle and determine operator size
             fun   = @(x,mode) opConvolve2Periodic_intrnl(fKernel,m,n,cflag,x,mode);
             nRows = m*n;
             nCols = m*n;
          elseif reflexive
             % ========= Reflexive convolution =========
             
             % Zero pad if needed
             bigKernel = zeros(m, n);
             bigKernel(1:size(kernel,1),1:size(kernel,2)) = kernel;
             kernel = bigKernel;
 
             % Precompute kernel in frequency domain
             e1 = zeros(size(kernel)); e1(1,1) = 1;
             fKernel = dct2(dctshift(kernel,offset))./dct2(e1);
 
             % Create function handle and determine operator size
             fun   = @(x,mode) opConvolve2Reflexive_intrnl(fKernel,m,n,cflag,x,mode);
             nRows = m*n;
             nCols = m*n;
          else
             % ========= Regular or truncated convolution =========
 
             % Shift kernel and add internal padding
             kernel = [kernel(offset(1):end,:);zeros(m-1,k(2));kernel(1:offset(1)-1,:)];
             kernel = [kernel(:,offset(2):end),zeros(size(kernel,1),n-1),kernel(:,1:offset(2)-1)];
             if truncated
                idx1 = 1:m;
                idx2 = 1:n;
             else
                idx1 = [size(kernel,1)-(offset(1)-2):size(kernel,1), 1:(m+k(1)-offset(1))];
                idx2 = [size(kernel,2)-(offset(2)-2):size(kernel,2), 1:(n+k(2)-offset(2))];
             end

             % Precompute kernel in frequency domain
             fKernel = fft2(full(kernel));

             % Create function handle and determine operator size
             fun   = @(x,mode) opConvolve2_intrnl(fKernel,k,m,n,idx1,idx2,cflag,x,mode);
             nRows = length(idx1) * length(idx2);
             nCols = m*n;
          end

          % Construct operator
          op = op@opSpot('Convolve2', nRows, nCols);
          op.cflag     = cflag;
          op.funHandle = fun;
       end % Constructor

    end % Methods       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = op.funHandle(x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef

%=======================================================================

function y = opConvolve2_intrnl(fKernel,k,m,n,idx1,idx2,cflag,x,mode)
if mode == 1
   fx = fft2([full(reshape(x,m,n)), zeros(m,k(2)-1); ...
              zeros(k(1)-1,n), zeros(k(1)-1,k(2)-1)]);
   y = ifft2(fKernel.*fx);
   y = y(idx1,idx2);
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
else
   z = zeros(m+k(1)-1,n+k(2)-1);
   z(idx1,idx2) = full(reshape(x,length(idx1),length(idx2)));
   y = ifft2(conj(fKernel).*fft2(z));
   y = y(1:m,1:n);
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
end
end

%======================================================================

function y = opConvolve2Periodic_intrnl(fKernel,m,n,cflag,x,mode)
if mode == 1
   y = ifft2(fKernel.*fft2(full(reshape(x,m,n))));
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
else
   y = ifft2(conj(fKernel).*fft2(full(reshape(x,m,n))));
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
end
end

%======================================================================

function y = opConvolve2Reflexive_intrnl(fKernel,m,n,cflag,x,mode)
if mode == 1
   y = idct2(fKernel.*dct2(full(reshape(x,m,n))));
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
else
   y = idct2(conj(fKernel).*dct2(full(reshape(x,m,n))));
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
end
end

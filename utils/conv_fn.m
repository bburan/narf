function newsig = conv_fn(signal, dim, fn, nwin, novl)
% CONV_FN: A function-application-based convolution along 1 dimension.
%
%  NEWSIG = conv_fn(SIGNAL, DIM, FN, NWIN, NOV);
%
% INPUTS:
%      SIGNAL    A 1D, 2D, 3D, or 4D matrix.
%      DIM       The dimension along which to apply the convolution.
%      FN        The function to apply during the downsampling. 
%                It should accept a vector and return a scalar.
%      NWIN      The number of elements of SIGNAL to bundle together.
%                Choose 10 to 'resample' 50Hz signal and make it 5Hz.
%      NOVP      The number of elements before and after the current point
%                to to overlap. Set to 0 if you don't want overlapping. 
%                Setting it to NWIN/2 will double the size of vectors
%                passed to FN. (The length would be NWIN/2 + NWIN + NWIN/2)
% 
% OPERATION:
% Input SIGNAL is grouped into NWIN-element long chunks along dimension
% DIM. Then NOVL elements groups are padded to the front and back of each
% chunk. If no such elements exist (such as at the beginning and end of
% SIGNAL), they will not be included in the chunk, so chunk length varies.
% 
% By this point, the chunk should have at most NOVL + NWIN + NOVL elements.
% This chunk is then passed to FN, which should be a function which accepts
% a vector. Useful picks for FN are @sum, @mean, or @max, corresponding to
% the L1, L2, and Linfty norms, but almost any function that accepts
% vectors should work.
%
% FN is executed as many times as there are chunks, and the return values 
% of FN are bundled together the convolution is returned. 
%
% If you set N to 1, you can create a windowed FIR response that is
% as nonlinear or arbitrary as you want, which can be useful. If you set it
% to 10 and use FN=@max, you are finding the maximum over each window. 
%
% OUTPUTS:
%      NEWSIG    The transformed SIGNAL matrix. 
%
% EXAMPLE:
%      x = rand([3,1000,2]);             % Create a random matrix
%      y1 = conv_fn(x, 2, @sum, 10, 5);  % Sum together groups of ten elements plus five on each side
%      hist(y1(1,:,1));                  % If a gaussian appears, you just demonstrated the central limit theorem.
%      plot(y1(1,2:end,1), y1(1,1:end-1,1), 'k.');  % Any groupings that appear demonstrate that nonzero NOVL resulted in each element correlating somewhat with the previous element
%      corrcoef(y1(1,2:end,1), y1(1,1:end-1,1));    % A correlation coefficient analysis confirms this statement
%
% 2012/10/31. Happy Halloween. Ivar Thorson.

if (nargin < 5)
    error('conv_fn() needs 5 or more arguments to work properly.');
end

dims_old = size(signal);  % Size of the old signal matrix
n_dims = ndims(signal);   % Number of dimensions of the matrix
n_old = dims_old(dim);    % Num of samples along the convolution dimension 

if n_dims > 4
    error('conv_fn() only works for 4D or smaller matrices');
end
if nwin < 1 | nwin ~= floor(nwin)
    disp(nwin);
    error('Sampling window size must be integer and 1 or greater.');
end
if novl < 0 | novl ~= floor(novl)
    error('Sampling window overlap must be integer and 0 or greater.');
end
if dim < 1 | dim > n_dims
    error('Convolution dimension is not indexable in the signal.');
end

% Prepare the new signal matrix
dims_new = dims_old;
if isequal(floor(dims_old(dim)/nwin), dims_old(dim)/nwin)
    dims_new(dim) = floor(dims_old(dim) / nwin);
else
    dims_new(dim) = floor((dims_old(dim) / nwin)) + 1;  
end
newsig = zeros(dims_new);

% Depending on which dimension was selected, filter that dimension
% When all you have is hammer, everything looks like a nail. 
% When all you have are for loops, the if conditional and a matrix data type...*sigh*
for a = 1:dims_new(1)
    if n_dims >= 2
        for b = 1:dims_new(2)
            if n_dims >= 3
                for c = 1:dims_new(3)
                    if n_dims >= 4
                        % Compute most bins
                        for d = 1:dims_new(4) 
                            if dim == 1
                                newsig(a,b,c,d) = fn(signal(max(1,(a-1)*nwin+1-novl):min([a*nwin+novl, n_old]),b,c,d)); 
                            elseif dim == 2
                                newsig(a,b,c,d) = fn(signal(a,max(1,(b-1)*nwin+1-novl):min([b*nwin+novl, n_old]),c,d)); 
                            elseif dim == 3
                                newsig(a,b,c,d) = fn(signal(a,b,max(1,(c-1)*nwin+1-novl):min([c*nwin+novl, n_old]),d)); 
                            elseif dim == 4
                                newsig(a,b,c,d) = fn(signal(a,b,c,max(1,(d-1)*nwin+1-novl):min([d*nwin+novl, n_old])));
                            end
                        end
                    else
                        if dim == 1
                            newsig(a,b,c) = fn(signal(max(1,(a-1)*nwin+1-novl):min([a*nwin+novl, n_old]),b,c));
                        elseif dim == 2
                            newsig(a,b,c) = fn(signal(a,max(1,(b-1)*nwin+1-novl):min([b*nwin+novl, n_old]),c));
                        elseif dim == 3 
                            newsig(a,b,c) = fn(signal(a,b,max(1,(c-1)*nwin+1-novl):min([c*nwin+novl, n_old])));
                        end
                    end
                end
            else
                if dim == 1
                    newsig(a,b) = fn(signal(max(1,(a-1)*nwin+1-novl):min([a*nwin+novl, n_old]), b));
                elseif dim == 2
                    newsig(a,b) = fn(signal(a,max(1,(b-1)*nwin+1-novl):min([b*nwin+novl, n_old])));
                end
            end
        end
    else
        newsig(a) = fn(signal(max(1,(a-1)*nwin+1-novl):min([a*nwin+novl, n_old])));
    end
end


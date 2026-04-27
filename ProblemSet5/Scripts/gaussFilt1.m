function x = gaussFilt1(x, sd, varargin)
% out = gaussFilt1(x, sd)
% out = gaussFilt1(x, sd, dim)
% 
% Filter 1 dimension of matrix or n-dimensional array "x" with a Gaussian
% of standard deviation sd. Default dimension to filter is dimension 1. To
% filter a different dimension, supply the dimension as argument "dim".
% Ends are treated symmetrically (see documentation for imfilter).
% 
% Matt Kaufman 2019-2020, University of Chicago

%% Optional argument

if isempty(varargin)
  dim = 1;
else
  dim = varargin{1};
end


%% Set up the filter

% Length of the filter to use. We'll make it 4 times the std dev, which is
% probably excessive
% Need to multiply by two because this is the full width, not the one-sided
% width
lenFactor = 4 * 2;
fLen = ceil(sd * lenFactor);

% Ensure the filter length is odd, which imgaussfilt needs to work
% correctly
fLen = fLen + mod(fLen + 1, 2);

% Set up the filter size for the first dimension and 1 for the second
% dimension (so we don't filter it!)
fSiz = [fLen 1];


%% Handle dim > 1
% If needed, permute x so that the first dimension is the dimension to
% filter. We'll just swap the dimension to filter with dimension 1, and
% swap back after. This is slightly inefficient if x is a matrix and dim ==
% 2, because imgaussfilt can handle this directly, but it's a small
% performance hit and made the next step much easier.

if dim > 1
  dimOrder = 1:ndims(x);
  dimOrder(1) = dim;
  dimOrder(dim) = 1;
  x = permute(x, dimOrder);
end


%% Handle ndims > 2
% Reshape into a matrix. This is ok, because we're always going to act on
% dimension 1 (if not, we've permuted so that we are)

siz = size(x);
x = reshape(x, siz(1), prod(siz(2:end)));


%% Do the actual filtering

x = imgaussfilt(x, sd, 'FilterSize', fSiz, 'Padding', 'symmetric');


%% Restore shape and dimension order

x = reshape(x, siz);

if dim > 1
  x = permute(x, dimOrder);
end

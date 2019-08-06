function y = nanstd(x,dim)
%NANSTD Standard deviation ignoring NaNs.
%   NANSTD(X) returns the same standard deviation treating NaNs 
%   as missing values.
%
%   For vectors, NANSTD(X) is the standard deviation of the
%   non-NaN elements in X.  For matrices, NANSTD(X) is a row
%   vector containing the standard deviation of each column,
%   ignoring NaNs.
%
%   NANSTD(X,DIM) takes the mean along the dimension DIM of X. 
%
%   See also NANMEAN, NANMEDIAN, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2000 The MathWorks, Inc. 
%   $Revision: 2.9 $  $Date: 2000/05/26 18:53:03 $

  
if isempty(x) % Check for empty input.
    y = NaN;
    return
end

if nargin == 1
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end


% Find mean
%avg = nanmean(x,dim);

tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

x = x - repmat(nanmean(x,dim),tile);  % Remove mean

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));
count = sum(~nans,dim);

% Protect against a scalar or column of all NaNs
i = find(count==0 | count == 1);
count(i) = ones(size(i)) + 1;
y = sqrt(sum(conj(x).*x,dim)./max(count-1,[],dim));
y(i) = i + NaN;

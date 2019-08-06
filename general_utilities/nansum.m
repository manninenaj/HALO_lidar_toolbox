function y = nansum(x,dim)
%NANSUM Sum ignoring NaNs.
%   NANSUM(X) returns the sum treating NaNs as missing values.  
%   For vectors, NANSUM(X) is the sum of the non-NaN elements in
%   X. For matrices, NANSUM(X) is a row vector containing the sum 
%   of the non-NaN elements in each column of X. 
%
%    See also NANMEDIAN, NANSTD, NANMIN, NANMAX, NANMEAN.

  
if nargin < 2; dim = 1; end
  
% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

% Protect against an entire column of NaNs
y = sum(x,dim);
i = find(all(nans,dim));
y(i) = i + NaN;

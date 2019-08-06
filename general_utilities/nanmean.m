function y = nanmean(x,dim)
%NANMEAN Average or mean ignoring NaNs.
%   NANMEAN(X) returns the average treating NaNs as missing values.  
%   For vectors, NANMEAN(X) is the mean value of the non-NaN
%   elements in X.  For matrices, NANMEAN(X) is a row vector
%   containing the mean value of each column, ignoring NaNs.
%
%   NANMEAN(X,DIM) takes the mean along the dimension DIM of X. 
%
%   See also NANMEDIAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2000 The MathWorks, Inc. 
%   $Revision: 2.10 $  $Date: 2000/05/26 18:53:02 $

if isempty(x) % Check for empty input.
    y = NaN;
    return
end

if nargin == 1
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end
  

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));
count = sum(~nans,dim);

% Protect against a column of all NaNs
i = find(count==0);
count(i) = ones(size(i));
y = sum(x,dim)./count;
y(i) = i + NaN;

function [result, count] = medianfilter(data, kernel, padvalue, options)
  
%   data = medianfilter(data, kernel, pad) 
%
% 2-D median filter with selectable kernel
% If function medfilt2 is available from the image processing 
% toolbox, will use that as it is faster (mex code) 
%
% Default kernel is [3 3]
%
% Default is no padding at the edges (nanmedian is used internally to
% give expected answer).
% If a value for padding is given, median is used internally.

if nargin < 1
  help medianfilter
  return;
end

if nargin < 2
  kernel = [3 3];
end

if nargin < 3
  padvalue = NaN;
  options = 'nopad';
else
  options = 'pad';
end



%%-- Check inputs --%%
if length(size(data)) > 2
  disp('This function has been written for 2-D matrices')
  return;
end


% check kernel...

%if prod(kernel)> 20
%  disp('may run out of memory..')
%end


 if exist('medfilt2')   
   result = medfilt2(data, kernel);
   count = [];
else


% Blindly assume that kernel should be symmetric about central point
% i.e. kernel should contain odd values.
kernel = (ceil(kernel ./ 2) .* 2) - 1;
padsize = floor(kernel ./2);
originalsize = size(data);

newsize = originalsize + 2.* padsize;
newmat = ones(newsize) .* padvalue;

% setup 3-D matrix according to kernel and padding 
largemat = repmat(newmat,[1 1 prod(kernel)]);
% fill in 3-D matrix with appropriately shifted 2-D data
for ii = 1:kernel(1)
  for jj = 1:kernel(2)
    tilenumber = (ii-1).*kernel(2) + jj;
    ix = ii:originalsize(1)+ii-1;
    iy = jj:originalsize(2)+jj-1;
    
    largemat(ix,iy,tilenumber) = data;
%  squeeze(largemat(:,:,tilenumber))
  end
end

selection = [padsize+1; padsize + originalsize];
ix = selection(1,1):selection(2,1);
iy = selection(1,2):selection(2,2);

switch options
 case 'nopad'
  result = nanmedian(largemat(ix,iy,:),3);
  count = nansum(largemat(ix,iy,:),3);
 otherwise
  result = median(largemat(ix,iy,:),3);
  count = sum(largemat(ix,iy,:),3);
end
end

function [Xo,nsamples,n_nonnan_samples] = windowSlider_v2(Xi,win,f,edges,prc_val,Yi)
%windowSlider_v2 can calculate mean, median, standard deviation, variance, min,
%max, skewness, kurtosis, fraction (inputs ones and zeros), or covariance 
%within a sliding window.
%
% Usage:
% Xo = windowSlider(Xi,win)
% Xo = windowSlider(Xi,win,f)
% Xo = windowSlider(Xi,win,f,edges)
% Xo = windowSlider(Xi,win,f,edges,prc_val) ... if 3rd input is @prctile
% Xo = windowSlider(Xi,win,f,edges,[],Yi) ... when calcuting covariance
% [Xo,nsamples] = windowSlider(Xi,...)
% [Xo,nsamples,n_nonnan_samples] = windowSlider(Xi,...)
%
% Inputs:
% - Xi                Input array, either a vector or a matrix
% - win               Size of sliding window, either scalar (3) or vector ([3,3])
% - f                 handle: @nanmean (default), @nanmedian, @nanvar, @nanstd, 
%                       @skewness, @kurtosis, @prctile, @fraction, @covariance.
% - edges             Edges padded with ('nans'), or with 'zeros'
% - prc_val           numeric scalar, specifies the percentile
% - Yi                2nd input array for covariance calculation, same dims as Xi
%
% Outputs:
% - Xo                Calculated quantity in the same dimensions as the input array
% - nsamples          Total number of samples
% - n_nonnan_samples  Number of non-nan samples
%             
%
% 12 October 2018
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% INAR, University of Helsinki, Finland

warning('Doensn''t work on large ''Xi'' and ''Yi'' arrays, runs out of memory!')

% Check the inputs and defaults
if nargin < 1
    error 'Not enough input arguments!'
elseif nargin < 2
    switch isvector(Xi)
        case 1
            win = 3;
        case 0
            win = [3,3];
    end
    f = @nanmean;
    edges = 'nans';
elseif nargin < 3
    f = @nanmean;
    edges = 'nans';
elseif nargin < 4
    edges = 'nans';
elseif nargin > 4 && isempty(edges)
    edges = 'nans';
elseif nargin > 6
    error 'Too many input arguments!'
end

% Check that the inputs meet the requirements
if ~isvector(Xi) && ~ismatrix(Xi)
    error('The first input has to be a 1-D or 2-D array!')
end
if ~isnumeric(win) && length(win) > 3 && ~(isscalar(win) || isvector(win))
    error(sprintf(['The second input specifies the length (if input array is a ' ...
        'vector) or size (if input array is a matrix) \nof the sliding ' ...
        'window, and has to be scalar or vector with max length of 2.']))
end
if length(win) > 1
    if abs(win(1)) > intmax
        error('The second input, size of the sliding window, has to be less than ''intmax''.')
    end
end
if ~isa(f, 'function_handle')
    error(sprintf(['The third input has to be a function handle: '...
        '@nanmean (default), @nanmedian, @nanvar, @nanstd,\n @skewness, @kurtosis, @prctile, @fraction, @covariance.']))
end

if ~ischar(edges) || ~(strcmp(edges,'nans') || strcmp(edges,'zeros'))
    error('The fourth input has to be a string specifying how the edges are handled: ''nans'' (default), or ''zeros''.')
end

if strcmp(func2str(f),'fraction')
    warning('With ''@fraction'' the input array has to be ones and zeros!')
end

% Determine padding
switch edges
    case 'nans'
        Xi_p = [nan(floor(win(1)/2),size(Xi,2)); Xi; nan(floor(win(1)/2),size(Xi,2))];
        Xi_p = [nan(size(Xi_p,1),floor(win(2)/2)), Xi_p, nan(size(Xi_p,1),floor(win(2)/2))];
        if nargin == 6
            Yi_p = [nan(floor(win(1)/2),size(Yi,2)); Yi; nan(floor(win(1)/2),size(Yi,2))];
            Yi_p = [nan(size(Yi_p,1),floor(win(2)/2)), Yi_p, nan(size(Yi_p,1),floor(win(2)/2))];
        end
    case 'zeros'
        Xi_p = [zeros(floor(win(1)/2),size(Xi,2)); Xi; zeros(floor(win(1)/2),size(Xi,2))];
        Xi_p = [zeros(size(Xi_p,1),floor(win(2)/2)), Xi_p, zeros(size(Xi_p,1),floor(win(2)/2))];
        if nargin == 6
            Yi_p = [zeros(floor(win(1)/2),size(Yi,2)); Yi; zeros(floor(win(1)/2),size(Yi,2))];
            Yi_p = [zeros(size(Yi_p,1),floor(win(2)/2)), Yi_p, zeros(size(Yi_p,1),floor(win(2)/2))];
        end
end

% Construct indices, which the sliding/running window would have
a0 = nan(1,size(Xi,2)); a0(1) = 0; for i = 2:size(Xi,2), a0(i) = a0(i-1)+size(Xi_p,1); end
rtmp1 = repmat(a0,1,size(Xi,1));
rtmp1s = sort(rtmp1);
a1 = repmat(rtmp1s,win(2)*win(1),1);
rtmp2 = repmat(1:size(Xi,1),1,size(Xi,2));
a2 = repmat(rtmp2,win(2)*win(1),1);
a3 = a1 + a2;
rtmp = repmat(transpose(0:win(1)-1),1,size(Xi,1)*size(Xi,2));
a4 = a3 + repmat(rtmp,win(2),1);
a5 = repmat(sort(repmat(transpose((1:win(2))-1)*size(Xi_p,1),win(1),1)),1,size(Xi,1)*size(Xi,2));
ind_running = a4 + a5;

% Calculate number of total samples and number of non-nan samples
n_nonnan_samples = reshape(sum(~isnan(Xi_p(ind_running))),size(Xi));
nsamples = reshape(repmat(size(ind_running,1),1,size(ind_running,2)),size(Xi));

switch func2str(f)
    case 'fraction'
    case 'prctile'
        if isempty(prc_val)
            error 'Specify percentile as the 5th input!'
        else
            Xo = reshape(prctile(Xi_p(ind_running),prc_val),size(Xi));
        end
    case 'skewness'
        Xo = reshape(my_skewness(Xi_p(ind_running)),size(Xi));
    case 'covariance'
        if nargin < 6
            error('For covariance the 6th input is required.')
        else
            Xo = reshape( nansum( ...
                ( ( Xi_p(ind_running) - repmat(nanmean(Xi_p(ind_running)),size(Xi_p(ind_running),1),1) ) .* ...
                  ( Yi_p(ind_running) - repmat(nanmean(Xi_p(ind_running)),size(Xi_p(ind_running),1),1) ) ) ) ./ ...
                repmat(size(ind_running,1),1,size(ind_running,2))  ,size(Xi) );
        end
    otherwise
        Xo = reshape(f(Xi_p(ind_running)),size(Xi));
end
end
% end

    function s = my_skewness(x,flag,dim)
        
        if nargin < 2 || isempty(flag)
            flag = 1;
        end
        if nargin < 3 || isempty(dim)
            % The output size for [] is a special case, handle it here.
            if isequal(x,[]), s = NaN; return; end
            
            % Figure out which dimension nanmean will work along.
            dim = find(size(x) ~= 1, 1);
            if isempty(dim), dim = 1; end
        end
        
        % Need to tile the output of nanmean to center X.
        tile = ones(1,max(ndims(x),dim));
        tile(dim) = size(x,dim);
        
        % Center X, compute its third and second moments, and compute the
        % uncorrected skewness.
        x0 = x - repmat(nanmean(x,dim), tile);
        s2 = nanmean(x0.^2,dim); % this is the biased variance estimator
        m3 = nanmean(x0.^3,dim);
        s = m3 ./ s2.^(1.5);
        
        % Bias correct the skewness.
        if flag == 0
            n = sum(~isnan(x),dim);
            n(n<3) = NaN; % bias correction is not defined for n < 3.
            s = s .* sqrt((n-1)./n) .* n./(n-2);
        end
    end

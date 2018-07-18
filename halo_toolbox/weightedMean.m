function y_wmean = weightedMean(y,y_error,method)
%weightedMean calculates mean weighted by the measurement errors, if 
%provided.
%
%Usage:
%y_wmean = weightedMean(y,y_error)
%y_wmean = weightedMean(y,y_error,method)  ---if method = 'unweighted' same as nanmean(y)
%y_wmean = weightedMean(y,[],method)       ---if method = 'unweighted' same as nanmean(y)
%y_wmean = weightedMean(y,[])              ---same as nanmean(y)
%
%Inputs:
%-y         numerical vector of values the mean is calculated from
%-error     numerical vector of uncertainties related to the values
%-method    string, either 'weighted' or 'unweighted', if 'y_error' is 
%           non-empty, but if empty, only 'unweighted' is accepted. default
%           'weighted'
%
%Outputs:
%-y_wmean    weighted mean    
%
%Created 2017-11-16
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

% Check inputs
if nargin < 2
    error('Not enough input arguments.')
end
if ~isnumeric(y(:)) && isscalar(y)
    error('1st input must be a numerical finite vector of matrix.')
end
if not(isempty(y_error))
    if ~isnumeric(y_error(:)) && isscalar(y_error) || ...
            size(y_error,1)~=size(y,1) || size(y_error,2)~=size(y,2)
        error(['2nd input must be a numerical finite vector of matrix,'...
            ' and has to have the same dimensions with the 1st input.'])
    end
end
if not(isempty(y_error)) && nargin < 3
    method = 'weighted'; % set default
elseif isempty(y_error) && nargin < 3
    method = 'unweighted';
elseif nargin == 3 && not(isempty(y_error))
    if not(ischar(method)) && not(any(strcmp({'weighted','unweighted'},method)))
        error('If provided, the 3rd input must a string and can be ''weighted'' or ''unweighted''.')
    end
elseif nargin == 3 && isempty(y_error) && strcmp(method,'weighted')
    error('Cannot calculate weighted mean if measurement error aren''t provided!')
elseif nargin > 3
    error('Too many inputs!')
end

if all(isnan(y(:)))
    y_wmean = nan;
else  
    switch method
        case 'weighted'
            % weights = uncertainties
            w = 1./y_error.^2;
            % V_p
            V_1 = nansum(w);% w^p, p = 1
        case 'unweighted' % for brevity; gives same result as 'mean(y)'
            nsamples = sum(~isnan(y));
            w = 1;
            V_1 = nsamples;
    end
    % quantity (iii)
    y_wmean = nansum(w .* y) ./ V_1;
end
end


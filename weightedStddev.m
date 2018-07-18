function y_stddev = weightedStddev(y,y_error,method)
%weightedStddev calculates either weighted or unweighted standard deviation
%unbiased by Gaussian noise and sample-size. Uses methods presented by Rimoldini (2014)
%http://dx.doi.org/10.1016/j.ascom.2014.02.001
%
%Usage:
%y_stddev = weightedStddev(y,y_error)
%y_stddev = weightedStddev(y,y_error,method)
%
%Inputs:
%-y         numerical vector of values the standard deviation is calculated from
%-y_error   numerical vector of uncertainties related to 'y'
%-method    string, either 'weighted' or 'unweighted', default ('weighted')
%
%Outputs:
%-y_stddev     standard deviation
%
%Created 2018-03-12
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

% Check inputs
if ~isnumeric(y(:)) || (~isvector(y) && ~ismatrix(y) && isscalar(y))
    error('1st input must be a numerical finite vector of matrix.')
end
if ~isnumeric(y_error(:)) || (~isvector(y_error) && ~ismatrix(y_error) && isscalar(y_error)) || ...
        size(y_error,1)~=size(y,1) || size(y_error,2)~=size(y,2)
    error(['2nd input must be a numerical finite vector of matrix,'...
        ' and has to have the same dimensions with the 1st input.'])
end
if nargin < 3
    method = 'weighted'; % set default
elseif nargin == 3
    if not(ischar(method)) && not(any(strcmp({'weighted','unweighted'},method)))
        error('If provided, the 3rd input must a string and can be ''weighted'' or ''unweighted''.')
    end
elseif nargin > 3
    error('Too many inputs!')
end

% NOTE: quantity and equation numbers refer to the Rimoldini (2014).

if all(isnan(y(:)))
    y_stddev = nan;
else    
    % Number of samples
    nsamples = sum(~isnan(y));
    nsamples(nsamples==0) = nan;
    switch method
        case 'unweighted'
            y_mean = nanmean(y);
            y_mean = repmat(y_mean,size(y,1),1);
            % quantity (iv), w = 1 and V_p = nsamples
            m_2 = nansum((y-y_mean).^2)./ nsamples;
            % Eq. (9), unweighted
            M_2 = nsamples ./ (nsamples-1) .* m_2;
            % sqrt of Eq.(28)
            y_stddev = real(sqrt(M_2 - 1./nsamples .* nansum(y_error.^2)));
        case 'weighted'
            % weights: uncertainties
            y_wmean = weightedMean(y,y_error,'weighted');
            y_wmean = repmat(y_wmean,size(y,1),1);
            w = 1./y_error.^2;
            % quantity (ii), V_p
            V_1 = nansum(w);% w^p, p = 1
            % weighted standard deviation, sqrt of Eq. (32)
            y_stddev = real(sqrt(nansum(w.*(y-y_wmean).^2./repmat(V_1,size(y,1),1)) - (nsamples-1)./V_1));
    end
end
end

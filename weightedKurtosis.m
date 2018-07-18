function y_kurto = weightedKurtosis(y,y_error,method)
%weightedKurtosis calculates either weighted or unweighted kurtosis unbiased by
%Gaussian noise and sample-size. Uses methods presented by Rimoldini (2014)
%http://dx.doi.org/10.1016/j.ascom.2014.02.001
%
%Usage:
%y_kurto = weightedKurtosis(y,y_error)
%y_kurto = weightedKurtosis(y,y_error,method)
%
%Inputs:
%-y         numerical vector of values the kurtosis is calculated from
%-y_error   numerical vector of uncertainties related to 'y'
%-method    string, either 'weighted' or 'unweighted', default ('weighted')
%
%Outputs:
%-y_kurto   kurtosis
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

if all(isnan(y(:)))
    y_kurto = nan;
else
    % Number of samples
    nsamples = sum(~isnan(y));
    nsamples(nsamples==0) = nan;
    switch method
        case 'unweighted'
            y_mean = nanmean(y);
            y_mean = repmat(y_mean,size(y,1),1);
            % quantity (iv)
            m_2 = nansum((y-y_mean).^2) ./ nsamples;
            m_4 = nansum((y-y_mean).^4) ./ nsamples;
            % Eq. (19)
            m_2_star = m_2 - ((nsamples - 1) ./ nsamples.^2) .* nansum(y_error.^2);
            % Eq. (21)
            m_4_star = m_4 - (6 .* (nsamples-2)) ./ nsamples.^2 .* ...
                nansum(y_error.^2 .* (y - y_mean).^2) - ...
                ((6 .* m_2_star) ./ (nsamples.^2)) .* nansum(y_error.^2) + ...
                (3 .* (nsamples-2).^2) ./ (nsamples.^3) .* nansum(y_error.^4) - ...
                (3 ./ nsamples.^4) .* (nansum(y_error.^2)).^2;
            % Eq. (22)
            m_2_2_star = (m_2_star.^2) - (4 ./ nsamples.^2) .* ...
                nansum(y_error.^2 .* (y - y_mean).^2) + ...
                (2 .* (nsamples-2)) ./ (nsamples.^3) .* nansum(y_error.^4) + ...
                (2 ./ nsamples.^4) .* (nansum(y_error.^2)).^2;
            % Eq. (9)
            M_2 = nsamples ./ (nsamples-1) .* m_2;
            % Eq.(28)
            K_2_star = M_2 - 1 ./ nsamples .* nansum(y_error.^2);
            % Eq. (31)
            K_4_star = (nsamples.^2 .* (nsamples + 1)) ./ ...
                ((nsamples-1) .* (nsamples-2) .* (nsamples-3)) .* m_4_star - ...
                (3 .* nsamples.^2) ./ ((nsamples-2) .* (nsamples-3)) .* m_2_2_star;
        case 'weighted'
            % weights: uncertainties
            y_wmean = weightedMean(y,y_error,'weighted');
            w = 1./y_error.^2;
            % quanttity (ii), V_p
            V_1 = nansum(w.^1);% p = 1
            V_2 = nansum(w.^2);% p = 2
            V_3 = nansum(w.^3);% p = 3
            V_4 = nansum(w.^4);% p = 4
            % calculate central moments when w = 1 / y_error^2
            % quantity (iv). m_r
            ymyw = nan(size(y)); for i = 1:size(y,1), ymyw(i,:) = y(i,:)-y_wmean; end
            m_2 = nansum(w.*(ymyw).^2) ./ V_1;% r = 2
            m_4 = nansum(w.*(ymyw).^4) ./ V_1;% r = 4
            % Eq. (32), '*' in the paper --> '_star' in here
            m_2_star = m_2 - (nsamples-1) ./ V_1;
            % Eq. (34)
            m_4_star = m_4 - 6./V_1 .* nansum((ymyw).^2 - y_error.^2./2) + ...
                (6.*m_2_star)./V_1 - 3./V_1.^2;
            % Eq. (35)
            m_2_2_star = m_2_star.^2 - (2.*(m_2_star+m_2)) ./ V_1;
            % Unbiased by sample-size and noise
            % Eq. (24)
            K_2_star = V_1.^2 ./ (V_1.^2-V_2) .* m_2_star;
            % Eq. (27)
            K_4_star = V_1.^2.*(V_1.^4-4.*V_1.*V_3+3.*V_2.^2) ./ ...
                ((V_1.^2-V_2).*(V_1.^4-6.*V_1.^2.*V_2+8.*V_1.*V_3+3.*V_2.^2-6.*V_4)) .* ...
                m_4_star - 3.*V_1.^2.*(V_1.^4-2.*V_1.^2.*V_2+4.*V_1.*V_3-3.*V_2.^2) ./ ...
                ((V_1.^2-V_2).*(V_1.^4-6.*V_1.^2.*V_2+8.*V_1.*V_3+3.*V_2.^2-6.*V_4)) .* m_2_2_star;
    end
    % Eq. (41)
    y_kurto = real(K_4_star./K_2_star.^2);
end
end

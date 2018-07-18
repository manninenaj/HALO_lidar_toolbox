function [i_outliers, i_outliers_high, i_outliers_low] =...
    findOutliers(x_outlr, y_outlr, adj)
%OUTLIERS function first carries out a robust linear regression to
% the input data. Then it calculates the leverage points and the
% indeces of the outliers based on the Cooks distance with
% optional adjustment provided by the user.
%
% Inputs:
% - x               x-values (vector), can contain nans
% - y               y-values (vector), can contain nans
% - adj             adjustment factor for Cooks distance threshold
%                   (scalar)
%
% Outputs:
% - iOutliers       indices of the outliers
% - iOutliersHigh   indices of the outliers, > fitted y-values
% - iOutliersLow    indices of the outliers, < fitted y-values

% Check number of inputs
if nargin ~= 3
    error 'Check number of inputs!'
end

% Check input validity
if isscalar(adj) ~= 1
    error 'Cook''s distance input has to be scalar!'
end

% Convert to column vectors
y_outlr = y_outlr(:); x_outlr = x_outlr(:);

% Select non nan from the original data
i_sel_org_outlr = find(~isnan(y_outlr)); % to be used later on
y_val_org_outlr = y_outlr(i_sel_org_outlr);
x_val_org_outlr = x_outlr(i_sel_org_outlr);

% Calculate the coefficients
b_outlr = my_robustfit(x_val_org_outlr(:),y_val_org_outlr(:));
[~, mID] = lastwarn; % iteration limit warning turned off
if ~isempty(mID)
    warning('off',mID)
end
p_coeff_outlr = [b_outlr(2) b_outlr(1)];

% Evaluate polynomial along the whole range
y_fit_outlr = polyval(p_coeff_outlr,x_outlr);

% Model matrix
X_outlr = [ones(length(x_val_org_outlr(:)),1), x_val_org_outlr(:)];

% Hat matrix
HatMatrix = X_outlr * ((X_outlr' * X_outlr) \ X_outlr');

% Leverage points
Leverage_points = diag(HatMatrix);

% Residuals
y_fit_valid_outlr = y_fit_outlr(~isnan(y_outlr)); % select non nan
y_resid_outlr = y_val_org_outlr(:) - y_fit_valid_outlr(:);

% Mean square error (MSE) of the predictor
MSE_outlr = nan(size(y_val_org_outlr));
for iMSE_outlr = 1:length(y_val_org_outlr)
    MSE_outlr(iMSE_outlr) = (y_fit_valid_outlr(iMSE_outlr) -...
        y_val_org_outlr(iMSE_outlr))^2;
end
MSE_outlr = sum(MSE_outlr) / numel(y_val_org_outlr);

% Cooks distance
D_Cook = nan(size(y_val_org_outlr));
for iCook = 1:length(y_val_org_outlr)
    D_Cook(iCook) = y_resid_outlr(iCook)^2 / ...
        (numel(p_coeff_outlr) * MSE_outlr) * ...
        (Leverage_points(iCook)/(1 - Leverage_points(iCook))^2);
end

% Initialize
i_outliers = i_sel_org_outlr;

% Indices where Cooks distance is less than threshold,
% here only the outliers remain while other are deleted
i_outliers(D_Cook < (4/numel(y_val_org_outlr)) * adj) = [];

% Indices where Cooks distance is less than threshold, and the
% observed values are larger than the fitted values
i_outliers_high = i_outliers(y_outlr(i_outliers) > ...
    y_fit_outlr(i_outliers));

% Indices where Cooks distance is less than threshold, and the
% observed values are lower than the fitted values
i_outliers_low  = i_outliers(y_outlr(i_outliers) < ...
    y_fit_outlr(i_outliers));

end


function [cloud_mask,signal_cld_scrd_outlr] = atmHALOsignalMasking(signal,range,win_size)
%atmoHALOsignalMasking creates a mask that can be used to filter
%atmoshpheric signal from Halo lidar measurements. Also, outputted are
%signal with everything but background filtered, filtered signal
%representing the full background field, and flags.
%
%
%
%
%
%


%-- Initial cloud-aerosol masking with 2d variance --%
% Divides the region of the calculated variance array contained in the
% furthest range bins (furthest 20%) into 64 subsections, which are used to
% find a dynamic threshold for the variance-based cloud-aerosol screening.

% Pad with nans to avoid the border effect
signal_pad_x = horzcat(signal, nan(round((win_size(1)-1)/2),...
    size(signal,1))');
signal_pad = horzcat(nan(round((win_size(1)-1)/2),...
    size(signal_pad_x,1))', signal_pad_x);

% 2D running variance with a window having dimensions equal to
% 'win_size'
signal_var_pad = nan(size(signal_pad));
for iC = 1+round((win_size(1)-1)/2):...
        size(signal_pad,2)-round((win_size(1)-1)/2)
    temp_array = signal_pad(:,iC-round((win_size(1)-1)/2):...
        iC+round((win_size(1)-1)/2));
    signal_var_pad(:,iC) = nanvar(temp_array,[],2);
end

% Remove padded nans
signal_var = signal_var_pad(:,round((win_size(1)-1)/2)+1:...
    end-round((win_size(1)-1)/2));
signal_var = real(log10(signal_var));
% Select variance in the furthest range bins (furthest 20%)
signal_var_up20 = signal_var(:,end-round(size(signal,2)*.2)+1:end);

% Find the nearest dimensions which are divisible by the number of
% provided subsections (default == 8^2 = 64) AND are larger than
% the original dimensions
nearest_div_mm = size(signal_var_up20,1);
while nearest_div_mm / 8 ~= ...
        fix(nearest_div_mm / 8)
    nearest_div_mm = nearest_div_mm + 1;
end
nearest_div_nn = size(signal_var_up20,2);
while nearest_div_nn / 8 ~= ...
        fix(nearest_div_nn / 8)
    nearest_div_nn = nearest_div_nn + 1;
end

% If needed, pad with zeros so that the variance array can be
% processed in equal-sized blocks
if (size(signal_var_up20,2)-nearest_div_nn) ~= 0
    signal_var_up20_pad_range = horzcat(signal_var_up20,...
        zeros(size(signal_var_up20,1),...
        nearest_div_nn-size(signal_var_up20,2)));
else
    signal_var_up20_pad_range = signal_var_up20;
end
if (size(signal_var_up20,1)-nearest_div_mm) ~= 0
    signal_var_up20_pad_all = vertcat(signal_var_up20_pad_range,...
        transpose(zeros(size(signal_var_up20_pad_range,2),...
        nearest_div_mm-size(signal_var_up20_pad_range,1))));
else
    signal_var_up20_pad_all = signal_var_up20_pad_range;
end

% Construct indeces for processing the array in blocks, such as:
%                 1 1 4 4 7 7
% block_indeces = 2 2 5 5 8 8
%                 3 3 6 6 9 9
ind_1st_col = sort(repmat(1:8,...
    1, size(signal_var_up20_pad_all, 2)/8));
ind_1st_block = repmat(ind_1st_col(:),...
    size(signal_var_up20_pad_all, 1)/8, 1);
ind_vec_all = ind_1st_block;
for i_sub_sect = 2:8
    ind_vec_all = vertcat(ind_vec_all, ind_1st_block + ...
        8 * (i_sub_sect - 1));
end
block_indeces = transpose(reshape(ind_vec_all,...
    [nearest_div_nn nearest_div_mm]));

% Remove the regions for padded zeros
if (abs(nearest_div_nn-size(signal_var_up20, 2))) ~= 0
    block_indeces(:, end - (nearest_div_nn -...
        size(signal_var_up20, 2)) + 1:end) = [];
end
if (abs(nearest_div_mm-size(signal_var_up20, 1))) ~= 0
    block_indeces(end - (nearest_div_mm -...
        size(signal_var_up20, 1)) + 1:end ,:) = [];
end

% Calculate mean variance in each block
block_mean = nan(1,max(block_indeces(:)));
for i_med = 1:max(block_indeces(:))
    block_mean(i_med) = ...
        nanmean(signal_var_up20_pad_all(block_indeces == i_med));
end

% Find the half of the blocks which have lowest median variance,
% i.e. which are least influenced by clouds and aerosols, and
% collect the calculated variance from those blocks into a
% reference array for finding the dynamic threshold
[~,sorted_block_ind] = sort(block_mean);
icond_refvar = ismember(block_indeces, sorted_block_ind(1:floor(...
    8 * 8/2)));
signal_var_ref = signal_var_up20(icond_refvar);
%         signal_var_ref = nan(size(signal_var_up20));
%         signal_var_ref(icond_refvar) = signal_var_up20(icond_refvar);

% Find dynamic threshold based on the number of screened pixels in
% the lowest most median variance regions of the signal. The
% allowed number of screened pixels is 10% of the total amount of
% pixels in the reference region, variable 'ratio'
th = 50; % initialise
i = 1;
signal_var_th = signal_var_ref;
signal_var_th(signal_var_ref > prctile(signal_var(:),th)) = nan;
ratio = sum(isnan(signal_var_th(:))) / numel(signal_var_ref);
while ratio(i) > .05
    i = i + 1;
    % Increase threshold every iteration
    th(i) = th(i-1) + 1;
    signal_var_th = signal_var_ref;
    signal_var_th(signal_var_ref > ...
        prctile(signal_var(:),th(i))) = nan;
    ratio(i) = sum(isnan(signal_var_th(:))) / ...
        numel(signal_var_ref);
    if ratio(i) == ratio(i-1) || th(i)>98
        break
    end
end

% Screen clouds based on dynamic variance threshold
signal_cloud_scr_var = signal;
signal_cloud_scr_var(signal_var > ...
    prctile(signal_var(:),th(i))) = nan;

% Final cloud-aerosol screening and filling
% Final cloud screening handles each vertical profile individually.
% Each vertical profile is fitted with a 1st degree polynomial.
% Note that the shape of the background can follow the shape of
% either 1st or 2nd degree polynomial. However, the magnitude of
% the change of the background as a function of range is
% insignificant compared to the magnitude of the outliers caused
% by the remnant cloud and aerosol signal. Thus, here the profiles
% are fitted with 1st degree polynomials. The outliers are removed
% with the 'findOutliers' function, which finds leverage points as
% well as points of high influence, and finally the indeces of the
% outliers based on their Cooks distance. Then, the regions which
% were cloud-screened have to be filled for the wavelet
% decomposition. The profiles are fitted with either 1st or 2nd
% degree polynomials based on the goodness-of-fit indicator,
% root-mean-square error.

% Initialize
signal_cld_scrd_outlr = signal_cloud_scr_var;
% Use the original uncorrected signal, i.e. no ripple removal, to
% construct the filled signal, which is used only in the wavelet
% transformation. If background *.txt files are absent, this has no
% effect since both original signal and ripple removed signal are
% in that case the same. Ripple removal passes the original as
% output if background *.txt files are absent.
for i_prof = 1:size(signal,1)
    if sum(~isnan(signal_cld_scrd_outlr(i_prof,:))) > ...
            numel(signal_cld_scrd_outlr(i_prof,:))*.05
        % Find outlier indices
        [~, i_outlrs_hi, ~] = findOutliers(range, ...
            signal_cld_scrd_outlr(i_prof,:), 1);        
        % Remove outliers
        signal_cld_scrd_outlr(i_prof,i_outlrs_hi) = nan;
    else
        % Assign as nans
        signal_cld_scrd_outlr(i_prof,:) = nan;
    end
end
% Cloud mask
cloud_mask = isnan(signal_cld_scrd_outlr);
end

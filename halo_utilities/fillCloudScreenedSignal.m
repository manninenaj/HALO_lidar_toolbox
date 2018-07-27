function [signal_fill,flag_nofit] = fillCloudScreenedSignal(signal,range,ignore_gates)
%fillCloudScreenedSignal The regions which were cloud-screened have to be 
%filled for the wavelet decomposition. The profiles are fitted with either 
%1st or 2nd degree polynomials based on the goodness-of-fit indicator,
%root-mean-square error.
%
%
%
%
%
%
%

signal_fill = nan(size(signal));
flag_nofit = nan(size(signal,1),1);
for i_prof = 1:size(signal,1)
    if sum(isnan(signal(i_prof,:))) ~= length(signal(i_prof,:))
        
        y_outlr_rmvd = signal(i_prof,:);
        x_outlr_rmvd = range;
        
        % Ignore additional number of range bins in case some
        % remnant aerosol signal did remain in the lower most range
        % bins. Note that the lower the fit reaches the more
        % accurate it is, tentatively.
        y_outlr_rmvd(1:ignore_gates) = nan;
        
        % Select only non nans
        y_outlr_rmvd_valid = y_outlr_rmvd(~isnan(y_outlr_rmvd));
        x_outlr_rmvd_valid = x_outlr_rmvd(~isnan(y_outlr_rmvd));
        
        % Get 1st and 2nd deg polynomial fits
        [B_1deg_final_scrn,stats_final_1deg] = ...
            my_robustfit(x_outlr_rmvd(:), y_outlr_rmvd(:));
        p_1deg_final_scrn = ...
            [B_1deg_final_scrn(2) B_1deg_final_scrn(1)];
        [B_2deg_final_scrn,stats_final_2deg] = ...
            my_robustfit([x_outlr_rmvd_valid(:) x_outlr_rmvd_valid(:).^2], y_outlr_rmvd_valid(:));
        p_2deg_final_scrn = ...
            [B_2deg_final_scrn(3) B_2deg_final_scrn(2) B_2deg_final_scrn(1)];
        
        % Evaluate for the whole range
        y_fit_final_scrn_1deg = polyval(p_1deg_final_scrn, range);
        y_fit_final_scrn_2deg = polyval(p_2deg_final_scrn, range);
        
        % Root Mean Squared Error (RMSE)
        RMSE_final_scrn_1deg = stats_final_1deg.ols_s;
        RMSE_final_scrn_2deg = stats_final_2deg.ols_s;
        
        % Initialize
        signal_fill(i_prof,:) = signal(i_prof,:);
        
        % Ignore additional number of range bins in case some
        % remnant aerosol signal did remain in the lower most range
        % bins. Though, the lower the fit reaches the more accurate
        % it is given the cloud screening was succesful.
        signal_fill(i_prof,1:ignore_gates) = nan;
        
        % Determine if theres enough data for a fit, based on
        % number of remaining points and location of the remaining
        % points; if too few data points OR if the data points
        % arent distributed sparsely enough, set flag --> 'no fit'
        signal_sel = ...
            signal(i_prof, 1:length(range) * .5);
        range_sel = range(1:length(range) * .5);
        if sum(~isnan(signal_sel)) < 5 || mean(range_sel(~isnan(signal_sel))) > prctile(range,45)
            
            % Do not fit, filled signal remains with nans within
            signal_fill(i_prof, isnan(signal(i_prof,:))) = nan;
            
            % Set flag
            flag_nofit(i_prof) = 1;
        else
            
            % If 1st degree polynomial fit is better
            if RMSE_final_scrn_1deg/RMSE_final_scrn_2deg < 1.1
                % Fill with 1st deg fits
                signal_fill(i_prof, isnan(signal_fill(i_prof,:))) = y_fit_final_scrn_1deg(isnan(signal_fill(i_prof,:)));
                
                % If 2nd degree polynomial fit is better
            else
                % Fill with 2nd deg fits
                signal_fill(i_prof, isnan(signal_fill(i_prof,:))) = y_fit_final_scrn_2deg(isnan(signal_fill(i_prof,:)));
            end
        end
    end
end
end



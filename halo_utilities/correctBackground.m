function [signal_corr, step_locations, flags, background] = ...
    correctBackground(signal, signal_orig, range, time, varargin)
%CORRECTBACKGROUND function corrects the background signal of the HALO
%Doppler lidar instrument. The background is corrected for step-changes and
%for the shape of the background within the step-changes respectively.
%
% Inputs:
%
% - signal              m-by-n matrix, the uncorrected SNR, 'm' equals the
%                       number of vertical profiles and 'n' equals the
%                       number of range bins (lengths of 'time' and 'range'
%                       inputs, respectively)
%
% - signal_orig         m-by-n matrix, the uncorrected SNR, 'm' equals the
%                       number of vertical profiles and 'n' equals the
%                       number of range bins (lengths of 'time' and 'range'
%                       inputs, respectively). No ripple removal. If ripple
%                       removal cant be done, should be the same as
%                       'signal'.
%
% - range               row or column vector, distances of each range bin
%                       from the instrument, length equals the number of
%                       range bins
%
% - time                row or column vector, time stamps of the vertical
%                       profiles, length equals the number of profiles
%
%
% Outputs:
%
% - signal_corrected    m-b-n matrix, the corrected SNR, dimensions equal
%                       the input 'signal'
%
% - step_locations      vector, indices of the detected step-changes in the
%                       background signal
%
% - flags               vector, flags for the vertical profiles
%                       0 = constrained fit through 1 (SNR)
%                       1 = 1st deg poly fit
%                       2 = 2nd deg poly fit
%                       3 = correction failed, replace with nans
%                       4 = correction failed, replaced with original
%
% - background          n-by-m matrix, the extracted uncorrected background
%                       signal, dimension equal the input 'signal'
%
% Additional arguments follow in the form of property/value pairs.
% Valid properties are: 'win_size', 'n_sub_sect', 'wavelet_level',
% 'correct_remnant'
%
% - 'win_size'          size of the sliding window, which is used in
%                       calculating the 2D variance in variance based cloud
%                       screening. DEFAULT: ('win_size', [33 1])
%
% - 'wavelet_level'     number of iterations for the wavelet decomposition.
%                       DEFAULT: ('wavelet_level',5)
%
% - 'correct_remnant'   determines how to handle the remnant outlier
%                       profiles, which are not presented well by the
%                       averaged  approach in section 3. DEFAULT:
%                       ('correct_remnant', 'original')
%                       'correct'   correct all of the profiles using
%                                   robust linear regression including the
%                                   remnant profiles
%                       'remove'    remove only the remnant outlier
%                                   profiles, replace with nans
%                       'original'  replace the remnant outlier profiles
%                                   with original profiles
%                       'none'      no correction applied
%
% - 'ignore'            By default, ignore (3) of the lowest range gates
%                       due to incontamination by the emitted pulse.
%                       Additionally, with the 'ignore' parameter user can
%                       assign e.g. 30 of the lowest most range gates to
%                       nans when calculating the shape of the background.
%                       Despite the cloud screening, some remnant aerosol
%                       signal might cause errors in calculating the shape.
%
% - 'cloud_mask'        cloud mask can be provided as an input as well. If
%                       not, then it is generated.
%
%
% Algorithm workflow (see more in detailed descriptions below)
% 0. Prepare data
% 1. Cloud screening
%      1.1 Crude cloud screening
%      1.2 Sub 2 km cloud screening
%      1.3 Final cloud screening
%      1.4 Calculation of the background shape
%      1.5 Filling the cloud-screened regions
% 2. Step detection
%      2.1 Multi-level 1D stationary wavelet decomposition
%      2.2 Peak detection i.e. locating step-changes
% 3. Correction for the step-changes and the shape of the background
% 4. OPTIONAL: Removal or correction of the remnant outlier profiles
%
% version 0.9.8
% 28 November 2016
% Antti Manninen
% antti.j.manninen@helsinki.fi

%% SET DEFAULTS
parameters.win_size        = [33 1];
parameters.wavelet_level   = 5;
parameters.correct_remnant = 'original';
parameters.ignore          = 3;
parameters.sizes           = [length(time) length(range)];
parameters.cloud_mask      = [];

%% CHECK THE INPUTS

% Was the 'parameters' struct supplied?
if ~isempty(varargin)
    % Check for overrides of the defaults
    parameters = parsePropertyValuePairs(parameters, varargin);
end

% Check whether the 'parameters' values can be accepted
parameters = checkParameters(parameters);

% Check the dimensions of the first three inputs
if ~ismatrix(signal) && ~isvector(time) && ~isvector(range) &&...
        length(time)  ~= size(signal,1) &&...
        length(range) ~= size(signal,2)
    error 'Check the input dimensions!'
else
    
    %% PREPARE THE DATA
    signal_0 = signal;
    range = range(:);
    time = time(:);
    % By default, ignore three of the lowest range gates due to
    % incontamination by the emitted pulse.
%     signal(:,1:4) = nan;
    signal(signal == 0) = nan;
    
    fprintf('\nStarting HALO background correction. This might take a while.\n')

    %% SCREENING AND FILLING
    if isempty(parameters.cloud_mask)
        [cloud_mask,~] = atmHALOsignalMasking(signal,range,parameters.win_size);
        % Use the original signal with ripples UNremoved, for the step detection 
        signal_cld_scrd_outlr = signal_orig;
        signal_cld_scrd_outlr(cloud_mask) = nan;
    else
        % Use the original signal with ripples UNremoved, for the step detection 
        signal_cld_scrd_outlr = signal_orig;
        signal_cld_scrd_outlr(parameters.cloud_mask) = nan;
    end
    
    fprintf('HALO background correction: filling signal for wavelet transformation...')
    signal_filled = fillCloudScreenedSignal(signal_cld_scrd_outlr,range,parameters.ignore);
    fprintf('done.\n')    

    fprintf('HALO background correction: wavelet transformation...')    
    %% MULTI-LEVEL 1D STATIONARY WAVELET DECOMPOSITION
    % The step-changes in the cloud-screened and filled signal are
    % detected by utilizing the stationary 1D wavelet decomposition
    % method. The output of the wavelet decomposition, detail
    % coefficients, describe the step-changes in the signal as peaks.
    % The detail coefficients are summed over each profile respectively
    % in order to make the peaks more pronounced and make the peak
    % detection more robust.
    
    % Use the default value 5 (recommended), OR use the value supplied
    % (suggest manual check)
    w_level = parameters.wavelet_level;
    
    % Initialize;
    detail_coeff_all  = [];
    i_switch_end   = 1; % for the zero padding
    
    % A range bin at a time
    for i_bin = 1:size(signal_filled, 2) * .75
        
        % Initialize
        zero_padded = signal_filled(:, end - i_bin);
        
        % The length of the zeropadded array ('zero_padded') has to be
        % divisible by (wavelet decomposition level)^2
        while round(size(zero_padded, 1) / 2^w_level) ~=...
                size(zero_padded,1 ) / 2^w_level
            
            % Pad signal with zeros to adjust its length
            zero_padded = [0; zero_padded; 0];
            
            % Remove one zero from the end and from the beginning
            % turn-by-turn
            if bitget(i_switch_end,1) % switch: odd or even
                
                % Remove zero from the end
                zero_padded(end) = [];
            else
                
                % Remove zero from the beginning
                zero_padded(1) = [];
            end
            
            % Add for iteration
            i_switch_end = i_switch_end + 1;
        end
        
        % Mark the locations of the added zeros
        cond_zeros = zero_padded == 0;
        
        % Convert zeros to nans for the wavelet function
        zero_padded(zero_padded == 0) = nan;
        
        % Call 'mySWT' function to calculate the stationary wavelet
        % decomposition up to the level indicated by 'w_level',
        % using the haar function for the convolution
        [~, detail_coeff] = mySWT(zero_padded, w_level, 'haar');
        
        % Leave only the highest level coefficients
        detail_coeff = detail_coeff(end,:);
        
        % Transpose to correct orientation
        detail_coeff = transpose(detail_coeff);
        
        % Get rid of nans which are padded zeros
        detail_coeff(cond_zeros) = [];
        
        % Concatenate arrays
        detail_coeff_all = horzcat(detail_coeff_all, detail_coeff);
    end
    
    % Sum up the coefficients along the bins
    detail_coeff_all_sum = nansum(abs(detail_coeff_all),2);
    
    % Convert zeros to nans
    detail_coeff_all_sum(detail_coeff_all_sum == 0) = nan;
    fprintf('done.\n')
    
    %% BACKGROUND STEP-CHANGE DETECTION
    % The step-changes in the signal are detected from inspecting the
    % detail coefficients from the output of the wavelet decomposition;
    % changes appear as peaks in the detail coefficients.
    % 'peakDetection' function detects peaks in signal based on given
    % delta value. The function searches for local minima and local
    % maxima from the signal (see more in detail description
    % 'peakDetection' function help).

    fprintf('HALO background correction: detecting steps in the background...')    
    % The steps occur as peaks (or outliers) in the detail coeffs.,
    % so look for peaks in a subset of data. Define what is a "peak"
    % i.e. difference between a "valley" and a "peak".
    [max_peaks_final_level, ~] = peakDetection(detail_coeff_all_sum.^2,...
        prctile(detail_coeff_all_sum, 75));
    fprintf('done.\n')

    fprintf('HALO background correction: performing initial correction...')
    if ~isempty(max_peaks_final_level)
        % Shift the locations based on the haar wavelet level impulse
        % distances
        step_locations = max_peaks_final_level(:, 1) + ...
            ((2^parameters.wavelet_level)/2)-1;
        step_locations(step_locations >= length(time) - 10) = [];
        step_locations([false; diff(step_locations)<10]) = [];
        %% CORRECT THE STEP-CHANGES AND THE SHAPE OF THE BACKGROUND
        % Removes the steps-changes and corrects the shape of the
        % background within the steps respectively. The correction is
        % carried out by fitting a surface to the cloud screened signal
        % within each step. The shape of the fitted surface is determined
        % by goodness-of-fit indicator root-mean-square error.
        
        % Initialize
        i_start_stp = 1;
        signal_shape_corrtd = nan(size(signal_0));
        background = nan(size(signal_0));
        flags = zeros(size(signal_0,1),1);
        
        for i_bkg_step = 1:length(step_locations) + 1
          
            % Determine background step end limit
            if i_bkg_step ~= length(step_locations) + 1
                i_end_stp = step_locations(i_bkg_step);
            else
                i_end_stp = size(signal, 1);
            end
            
            %% DRIFT CORRECTION
            % Median signal in range bins per each time stamp within a step
            med_signal_stp = nanmedian(signal_cld_scrd_outlr(...
                i_start_stp:i_end_stp,:),2);
            
            % Exclude nans
            x_bkg_stp = 1:length(time(i_start_stp:i_end_stp));
            x_bkg_stp(isnan(med_signal_stp)) = [];
            med_signal_stp(isnan(med_signal_stp)) = [];
            
            % Calculate linear fit and evaluate for the whole step
            linear_coeffs = my_robustfit(x_bkg_stp(:), med_signal_stp(:));
            p_coeff_drift = [linear_coeffs(2) linear_coeffs(1)];
            step_fitted = polyval(p_coeff_drift, ...
                1:length(time(i_start_stp:i_end_stp)));
            
            % Correct the drift with a step
            signal_drift_corrtd = signal_0(i_start_stp:i_end_stp,:) -...
                repmat(step_fitted(:),1,length(range)) +...
                p_coeff_drift(end);
            
            %% CORRECT BACKGROUND SHAPE AND STEP CHANGES
            % Determines the shape of the background between the
            % step-changes, and corrects for the shape and for the steps in
            % the background signal.
            
            % Exclude clouds
            signal_drift_corrtd(isnan(signal_cld_scrd_outlr(...
                i_start_stp:i_end_stp,:))) = nan;
            
            % Median signal per range bin within a step
            signal_med_bin = nanmedian(signal_drift_corrtd);
            
            % Determine is there enough data for higher order fit, depends
            % on the number and location of the remaining points. If there
            % are too few data points OR if the data points arent
            % distributed sparsely enough, set flag --> no higher order
            % fit
            Signal_MedStep_ySel = signal_med_bin(1:length(range)*.5);
            range_cld_scrd_outlr_sel = range(1:length(range)*.5);
            if sum(~isnan(Signal_MedStep_ySel)) < 5 || ...
                    mean(range_cld_scrd_outlr_sel(...
                    ~isnan(Signal_MedStep_ySel))) > prctile(range,45)
                
                % Set flag
                flag_fit_step = 1;
            else
                flag_fit_step = 0;
            end
            
            % Select only non nans
            ysel = signal_drift_corrtd;
            xsel = transpose(repmat(range(:),1,size(ysel,1)));
            [xsel_sort,isort] = sort(xsel(:));
            ysel = ysel(:); ysel_sort = ysel(isort);
            y_final_valid = ysel_sort(~isnan(ysel_sort));
            x_final_valid = xsel_sort(~isnan(ysel_sort));
            
            % Calculate 1st, 2nd deg, and constrained polynomial fits
            [B_prof_1deg,stats_1deg] = my_robustfit(x_final_valid(:), ...
                y_final_valid(:));
            [B_prof_2deg,stats_2deg] = my_robustfit([x_final_valid(:) ...
                x_final_valid(:).^2], y_final_valid(:));
            p_1deg_prof = [B_prof_1deg(2) B_prof_1deg(1)];
            p_2deg_prof = [B_prof_2deg(3) B_prof_2deg(2) B_prof_2deg(1)];
            % Evaluate coefficients
            y_fit_prof_1deg = polyval(p_1deg_prof, range);
            y_fit_prof_2deg = polyval(p_2deg_prof, range);

            % If 1st degree polynomial fit is better
            if stats_1deg.ols_s/stats_2deg.ols_s < 1.1 || ...
                    flag_fit_step == 1
                
                % Correct for step change and background shape
                signal_shape_corrtd(i_start_stp:i_end_stp,:) = ...
	            signal_drift_corrtd - repmat(transpose(y_fit_prof_1deg(:)),...
                    length(i_start_stp:i_end_stp), 1) + ...
                    p_1deg_prof(2);
                
                % Set flag to '1' if 1st deg fit was used
                flags(i_start_stp:i_end_stp,1) = ...
                    ones(length(i_start_stp:i_end_stp), 1);
                
                % Combine corrections to form corrected background
                background(i_start_stp:i_end_stp,:) = ...
                    (repmat(step_fitted(:),1,length(range)) - ...
                    p_coeff_drift(end))...
	      + repmat(transpose(y_fit_prof_1deg(:)), ...
                    length(i_start_stp:i_end_stp), 1);
                
                % If 2nd degree polynomial fit is better
            else
                % Correct for step change and background shape
                signal_shape_corrtd(i_start_stp:i_end_stp,:) = ...
                    signal_0(i_start_stp:i_end_stp,:) -...
 	            repmat(transpose(y_fit_prof_2deg(:)), ...
                    length(i_start_stp:i_end_stp), 1) + ...
                    p_2deg_prof(3);
                
                % Set flag to '2' 2nd deg fit was used
                flags(i_start_stp:i_end_stp,1) = ...
                    repmat(2,length(i_start_stp:i_end_stp), 1);
                
                % Combine corrections to form corrected background
                background(i_start_stp:i_end_stp,:) = ...
                    (repmat(step_fitted(:),1,length(range)) - ...
                    p_coeff_drift(end))...
	            + repmat(transpose(y_fit_prof_2deg(:)), ...
                    length(i_start_stp:i_end_stp), 1);
            end
            
            % For iteration
            i_start_stp = i_end_stp + 1;
        end
        
        % Correct the step changes and the shape of the background
        signal_shape_corrtd = signal_0 - background + 1;
        fprintf('done.\n')
    else
        % If the step detection failed, skip it
        signal_shape_corrtd = signal_0;
        step_locations = [];
        flags = repmat(4,size(signal_0));
        fprintf('failed.\n')
    end
    %% REMOVAL OR CORRECTION OF REMNANT OUTIER PROFILES
    % For instruments that are not operating optimally, the background
    % noise for some profiles may not be very well represented by the
    % the averaged approach. The outlier profiles can then be flagged
    % and rejected, or the user may choose to apply the background
    % noise profile shape detection and correction on a
    % profile-by-profile basis.
    fprintf('HALO background correction: performing final correction...')
    

    % Initialise
    signal_remn = nan(size(signal_0));
    
    switch parameters.correct_remnant
        % 'correct': correct all of the profiles using robust linear
        %            regression including the remnant profiles
        % 'remove':  remove only the remnant outlier profiles
        % 'original':    first correct all profiles using robust linear
        %            regression, and then remove the outlier profiles
        % 'none':    no correction
        case 'correct'
            fprintf('profile-by-profile...')
            for i_remn = 1:size(signal_0,1)
                % Select only background signal and apply cloudmask
                % calculated from the shape corrected signal
                y_remn = signal_shape_corrtd(i_remn,:);
                y_remn(parameters.cloud_mask(i_remn,:)) = nan;
                y_remn(1:parameters.ignore) = nan;
                if sum(isnan(y_remn)) ~= length(y_remn)
                    y_r_val = y_remn(~isnan(y_remn));
                    x_r_val = range(~isnan(y_remn));
                    
                    % Calculate robust bisquare linear fit
                    b_remn = my_robustfit(x_r_val(:),y_r_val(:));
                    if isempty(b_remn)
                        % No fit, no correction
                        signal_remn(i_remn,:) = signal_shape_corrtd(i_remn,:);
                    else
                        p_c_remn = [b_remn(2) b_remn(1)];
                        y_f_r = polyval(p_c_remn,range);
                        signal_remn(i_remn,:) = signal_shape_corrtd(i_remn,:)-transpose(y_f_r(:))+1;
                    end
                end
            end
            
            % Final corrected signal
            signal_corr = signal_remn;

        case 'remove'
            fprintf('by replacing "bad" profiles with NaNs...')
            
            for i_remn = 1:size(signal_0,1)
                % Select only background signal and apply cloudmask
                % calculated from the shape corrected signal
                y_remn = signal_shape_corrtd(i_remn,:);
                y_remn(cloud_mask(i_remn,:)) = nan;
                y_remn(1:parameters.ignore) = nan;
                if sum(isnan(y_remn)) ~= length(y_remn)
                    y_r_val = y_remn(~isnan(y_remn));
                    x_r_val = range(~isnan(y_remn));
                    
                    % Calculate robust bisquare linear fit
                    b_remn = my_robustfit(x_r_val(:),y_r_val(:));
                    p_c_remn = [b_remn(2) b_remn(1)];
                    y_f_r = polyval(p_c_remn,range);
                    signal_remn(i_remn,:) = ...
                        signal_shape_corrtd(i_remn,:) - transpose(y_f_r(:)) + 1;
                end
            end
            
            % Initialize with corrected signal
            signal_remn_cld_scrd = signal_remn;
            
            % Remove clouds, take median for each profile
            signal_remn_cld_scrd(isnan(signal_cld_scrd_outlr)) = nan;
            signal_remn_cld_scrd(:,1:parameters.ignore) = nan;
            signal_remn_peaks = nanmedian(signal_remn_cld_scrd,2);
            
            % Find the remant profiles with peak detection from the
            % median signal in profiles by using a 6 x std dev as peak
            % threshold
            [remn_peaks_max, remn_peaks_min] = ...
                peakDetection(signal_remn_peaks, ...
                6 * nanstd(signal_remn_peaks));
            
            % Collect min and max peak locations
            if ~isempty(remn_peaks_max) && ~isempty(remn_peaks_min)
                remn_peak_loc = vertcat(remn_peaks_max(:,1), ...
                    remn_peaks_min(:,1));
                
                % Set outlier profiles to nan
                signal_remn(remn_peak_loc,:) = nan;
                flags(remn_peak_loc) = 3;
            end
            
            % Final corrected signal
            signal_corr = signal_remn;

        case 'original'
            fprintf('by leaving "bad" profiles untouched...')
            for i_remn = 1:size(signal_0,1)             
                % Select only background signal and apply cloudmask
                % calculated from the shape corrected signal
                y_remn = signal_shape_corrtd(i_remn,:);
                y_remn(cloud_mask(i_remn,:)) = nan;
                y_remn(1:parameters.ignore) = nan;
                if sum(isnan(y_remn)) ~= length(y_remn)
                    y_r_val = y_remn(~isnan(y_remn));
                    x_r_val = range(~isnan(y_remn));
                    
                    % Calculate robust bisquare linear fit
                    b_remn = my_robustfit(x_r_val(:),y_r_val(:));
                    [~, mID_r] = lastwarn; % iteration limit warning off
                    if ~isempty(mID_r), warning('off',mID_r), end
                    p_c_remn = [b_remn(2) b_remn(1)];
                    y_f_r = polyval(p_c_remn,range);
                    signal_remn(i_remn,:) = ...
                        signal_shape_corrtd(i_remn,:) - transpose(y_f_r(:)) + 1;
                end
            end
            
            % Initialize with corrected signal
            signal_remn_cld_scrd = signal_remn;
            
            % Remove clouds, take median for each profile
            signal_remn_cld_scrd(isnan(signal_cld_scrd_outlr)) = nan;
            signal_remn_cld_scrd(:,1:parameters.ignore) = nan;
            signal_remn_peaks = nanmedian(signal_remn_cld_scrd,2);
            
            % Find the remant profiles with peak detection from the
            % median signal in profiles by using a 6 x std dev as peak
            % threshold
            [remn_peaks_max, remn_peaks_min] = ...
                peakDetection(signal_remn_peaks, ...
                6 * nanstd(signal_remn_peaks));
            
            % Collect min and max peak locations
            if ~isempty(remn_peaks_max) && ~isempty(remn_peaks_min)
                remn_peak_loc = vertcat(remn_peaks_max(:,1), ...
                    remn_peaks_min(:,1));
                
                % Set outlier profiles to original
                signal_remn(remn_peak_loc,:) = signal(remn_peak_loc,:);
                flags(remn_peak_loc) = 4;
            end
            
            % Final corrected signal
            signal_corr = signal_remn;
            
        otherwise
            % includes option none %
            signal_corr = signal_shape_corrtd;
    end
    background = signal_0 - signal_corr;
    fprintf('done.\n')
end

% Subfunction - parsePropertyValuePairs
%DErrico, John (2006). Parsing property/value pairs for function input
%(http://www.mathworks.com/matlabcentral/fileexchange/9082-parse-pv-pairs),
%MATLAB Central File Exchange. Retrieved Oct 7, 2015.

    function params = parsePropertyValuePairs(params, pv_pairs)
        %PARSEPROPERTYVALUEPAIRS: parses sets of property value pairs
        % usage: params = parse_pv_pairs(default_parameters, pv_pairs)
        %
        % arguments: (input)
        %  default_parameters - structure, with one field for every
        %   potential property/value pair. Each field will contain the
        %   default value for that property. If no default is supplied for
        %   a given property, then that field must be empty.
        %
        %  pv_array - cell array of property/value pairs.
        %   Case is ignored when comparing properties to the list of field
        %   names. Also, any unambiguous shortening of a field/property
        %   name is allowed.
        %
        % arguments: (output)
        %  params - parameter struct that reflects any updated
        %  property/value pairs in the pv_array.
        %
        % Example usage:
        % First, set default values for the parameters. Assume we
        % have four parameters that we wish to use optionally in
        % the function examplefun.
        %
        %  - 'viscosity', which will have a default value of 1
        %  - 'volume', which will default to 1
        %  - 'pie' - which will have default value 3.141592653589793
        %  - 'description' - a text field, left empty by default
        %
        % The first argument to examplefun is one which will always be
        % supplied.
        %
        %   function examplefun(dummyarg1,varargin)
        %   params.Viscosity = 1;
        %   params.Volume = 1;
        %   params.Pie = 3.141592653589793
        %
        %   params.Description = '';
        %   params=parse_pv_pairs(params,varargin);
        %   params
        %
        % Use examplefun, overriding the defaults for 'pie', 'viscosity'
        % and 'description'. The 'volume' parameter is left at its default.
        %
        %   examplefun(rand(10),'vis',10,'pie',3,'Description',
        %    'Hello world')
        %
        % params =
        %     Viscosity: 10
        %        Volume: 1
        %           Pie: 3
        %   Description: 'Hello world'
        %
        % Note that capitalization was ignored, the property 'viscosity'
        % was truncated as supplied. Also, note that the order the pairs
        % were supplied was arbitrary.
        
        npv = length(pv_pairs);
        n = npv/2;
        
        if n~=floor(n)
            error 'Property/value pairs must come in PAIRS.'
        end
        if n<=0
            % just return the defaults
            return
        end
        
        if ~isstruct(params)
            error 'No structure for defaults was supplied'
        end
        
        % there was at least one pv pair. process any supplied
        propnames = fieldnames(params);
        lpropnames = lower(propnames);
        for i_sf=1:n
            p_i = lower(pv_pairs{2*i_sf-1});
            v_i = pv_pairs{2*i_sf};
            
            ind = strmatch(p_i,lpropnames,'exact');
            if isempty(ind)
                ind = find(strncmp(p_i,lpropnames,length(p_i)));
                if isempty(ind)
                    error(['No matching property found for: ',...
                        pv_pairs{2*i_sf-1}])
                elseif length(ind)>1
                    error(['Ambiguous property name: ',pv_pairs{2*i_sf-1}])
                end
            end
            p_i = propnames{ind};
            
            % override the corresponding default in params
            params = setfield(params,p_i,v_i); %#ok
            
        end
    end

%% Subfunction - checkParameters
    function params = checkParameters(params)
        %CHECKPARAMETERS checks that the input parameters are correct
        
        % amount of range bins to be ignored at the bottom
        if isempty(params.ignore)
            params.ignore = 4;
        elseif not(isnumeric(params.ignore) && ...
                isscalar(params.ignore) && ...
                params.ignore <params.sizes(2) * .5)
            error(['''ignore'' parameter has to be numeric scalar and' ...
                ' cannot be larger than 1/2 times the number of' ...
                ' range gates. Default is %d'],3)
        end
        
        % win_size == 1 by default
        if isempty(params.win_size)
            params.win_size = [33 1];
        else
            if not(isnumeric(params.win_size) && ...
                    isvector(params.win_size) && ...
                    numel(params.win_size) == 2)
                error(['''win_size'' parameter must be a numeric' ...
                    ' vector with length of two, default [33,1]'])
            end
        end
        
        % cloud_mask = [] by default
        if isempty(params.cloud_mask)
            params.cloud_mask = [];
        else
            if not(islogical(params.cloud_mask)) && all(size(params.cloud_mask,1) == params.sizes)
                error(['''cloud_mask'' parameter must be a logical array' ...
                    ' with same dimensions as input signal.'])
            end
        end
        
        % correct_remnant - must be one of 3 options
        valid = {'correct', 'remove', 'original', 'none'};
        if isempty(params.correct_remnant)
            % default == 'original'
            params.correct_remnant = 'original';
        end
        ind = find(strncmpi(params.correct_remnant,valid,...
            length(params.correct_remnant)));
        if (length(ind) == 1)
            params.correct_remnant = valid{ind};
        else
            error(['Invalid input for remnant correction parameter.' ...
                ' Valid options are: ''%s'', ''%s'', ''%s''' ...
                '(default), or ''%s'''],valid{1},valid{2},valid{3},...
                valid{4})
        end
        
        % wavelet_level == 5 by default
        if isempty(params.wavelet_level)
            params.wavelet_level = 5;
        elseif not(isnumeric(params.wavelet_level) && ...
                isscalar(params.wavelet_level)) && ...
                (length(params.wavelet_level)>1)
            error ''wavelet_level' must be a numeric scalar. Default is 5.'
        end
        
    end


%-- Subfunction - mySWT
    function [approx_coeff, detail_coeff] = mySWT(x, w_level, wavelet_type)
        % Stationary wavelet decomposition for correcting background
        %
        % USAGE:
        %
        %   [approx_coeff, detail_coeff] = my_swt(x,w_level,wavelet_type);
        %
        % INPUTS
        %
        %    MY_SWT(X,LEVEL,'wavelet_type') computes the stationary wavelet
        %    decomposition of the signal X up to LEVEL, using wavelet of 
        %    type 'wavelet_type'. LEVEL must be a positive integer and X 
        %    should have a length divisible by 2^LEVEL.
        %
        %    Currently, 'Haar' is the only supported wavelet type.
        %
        %
        % OUTPUTS
        %
        %   Outputs are the approximation and detailed coefficients for 
        %   each level.
        %
        
        % Use row vector.
        x = transpose(x(:));
        lx = length(x);
        % Check that the length of x is divisible by 2^(wavelet 
        % decomposition level)
        if rem(lx,2^w_level)>0
            disp(['SWT input should have length of divisible by ' ...
                '2^(wavelet decomposition level).' ...
                ' Something wrong when zero padding']);
            return
        end
        
        % Get decomposition filters.
        switch wavelet_type
            case 'haar'
                lopass = [ 1./sqrt(2) 1./sqrt(2)];
                hipass = [-1./sqrt(2) 1./sqrt(2)];
            otherwise
                disp('Ask for new version');
                return
        end
        
        % Compute stationary wavelet coefficients.
        approx_coeff = zeros(w_level,lx);
        detail_coeff = zeros(w_level,lx);
        
        for k = 1:w_level
            
            % Extension
            lf = length(lopass); % length of filter
            
            % Discrete Wavelet Transform mode is periodisation
            x = extendPeriodDWT(x,lf/2);
            
            % Decomposition
            detail_coeff(k,:) = extractVector(conv2(x(:)',hipass(:)',...
                'full'),lx,lf+1);
            approx_coeff(k,:) = extractVector(conv2(x(:)',lopass(:)',...
                'full'),lx,lf+1);
            
            % Dyadic upsampling of filters
            tmp = zeros(1,2.*lf);
            tmp(1:2:2 * lf) = lopass;
            lopass = tmp;
            tmp = zeros(1,2.*lf);
            tmp(1:2:2 * lf) = hipass;
            hipass = tmp;
            
            % Update x
            x = approx_coeff(k,:);
        end
    end

%-- Subfunction - extractVector
    function y = extractVector(x, len, start)
        %EXTRACTVECTOR extracts a vector from within a larger vector
        
        y = x;
        finish = start + len - 1;
        y = y(start:finish);
        
    end

%-- Subfunction - extendPeriodDWT
    function x = extendPeriodDWT(x,lf)
        %EXTENDPERIODDWT extends the DWT using periodisation
        
        length_of_x = length(x);
        if rem(length_of_x,2)
            x(length_of_x+1) = x(length_of_x);
            length_of_x = length_of_x+1;
        end
        
        I = [length_of_x-lf+1:length_of_x , 1:length_of_x , 1:lf];
        if length_of_x<lf
            I = mod(I,length_of_x);
            I(I==0) = length_of_x;
        end
        
        x = x(I);
        
    end

%-- Subfunction - peakDetection
% Modified after Billauer, Eli (2012). Peak detection using MATLAB
% (http://www.billauer.co.il/peakdet.html). Last accessed 7 Oct 2015

    function [max_tab_peakd, min_tab_peakd] = peakDetection(vect_peakd,...
            delta_peakd, ex_peakd)
        %PEAKDETECTION Detects peaks in a vector
        % [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local maxima and 
        % minima ("peaks") in the vector V. MAXTAB and MINTAB consists of 
        % two columns. Column 1 contains indices in V, and column 2 the 
        % found values.
        %
        % With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices in 
        % MAXTAB and MINTAB are replaced with the corresponding X-values.
        %
        % A point is considered a maximum peak if it has the maximal value,
        % and was preceded (to the left) by a value lower by DELTA.
        
        max_tab_peakd = [];
        min_tab_peakd = [];
        
        vect_peakd = vect_peakd(:); % in case this wasnt a proper vector
        
        % Check inputs
        if nargin < 3
            ex_peakd = transpose(1:length(vect_peakd));
        else
            ex_peakd = ex_peakd(:);
            if length(vect_peakd)~= length(ex_peakd)
                error(['Input vectors ''vect_peakd'' and ''ex_peakd'' ' ...
                    'must have the same length']);
            end
        end
        
        if (length(delta_peakd(:)))>1
            error('Input argument ''delta_peakd'' must be a scalar');
        end
        
        if delta_peakd <= 0
            error('Input argument ''delta_peakd'' must be positive');
        end
        
        % Prepare data
        mn_peakd = Inf;
        mx_peakd = -Inf;
        mnpos_peakd = nan;
        mxpos_peakd = nan;
        
        lookformax_peakd = 1;
        
        for i_peakd = 1:length(vect_peakd)
            this_peakd = vect_peakd(i_peakd);
            if this_peakd > mx_peakd
                mx_peakd = this_peakd;
                mxpos_peakd = ex_peakd(i_peakd);
            end
            if this_peakd < mn_peakd
                mn_peakd = this_peakd;
                mnpos_peakd = ex_peakd(i_peakd);
            end
            
            if lookformax_peakd
                if this_peakd < mx_peakd-delta_peakd
                    max_tab_peakd = [max_tab_peakd; mxpos_peakd mx_peakd];
                    mn_peakd = this_peakd; mnpos_peakd = ex_peakd(i_peakd);
                    lookformax_peakd = 0;
                end
            else
                if this_peakd > mn_peakd+delta_peakd
                    min_tab_peakd = [min_tab_peakd; mnpos_peakd mn_peakd];
                    mx_peakd = this_peakd; mxpos_peakd = ex_peakd(i_peakd);
                    lookformax_peakd = 1;
                end
            end
        end
    end
end

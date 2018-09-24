function calculateHALOatmBLclassificationProduct(site,DATES,dt,dtskewn,weighting)
%calculateHALOatmBLclassificationProduct generates boundary layer classification product
%from calculated Doppler lidar quantities, and writes the product into a daily *.nc file.
%
% Usage:
% calculateHALOatmBLclassificationProduct(site,DATES)
% calculateHALOatmBLclassificationProduct(site,DATES,dt)
% calculateHALOatmBLclassificationProduct(site,DATES,dt,dtskewn)
% calculateHALOatmBLclassificationProduct(site,DATES,dt,dtskewn,weighting)
% calculateHALOatmBLclassificationProduct(site,DATES,[],[]  ,weighting)
%
% Inputs:
% -site          string, site name, e.g. site = 'kuopio'
% -DATES         scalar or vector, numeric, e.g. DATES = 20170401
%                or DATES = [20170401 20170431]
% -dt            scalar or vector, numeric, temporal resolution in minutes, 
%                e.g. dt = [1 3 5], or dt = 30, by default dt = [3 5 10 30]
% -dtskewn       scalar or vector, numeric, temporal resolution in minutes, 
%                e.g. dt = 20, by default dt = 30.
% -weighting     logical 'true' or 'false'
%
%
% Created 2018-03-08
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

% Check inputs
if nargin < 2
  error('''site'', and ''DATES'' are required inputs!')
end
if ~ischar(site)
  error('The first input ''site'' must be a string.')
end
if length(DATES)>2
    error('''DATES'' can have max. length of 2.')
elseif length(DATES)==1
    DATEstart = DATES; DATEend = DATES;
elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
       length(num2str(DATES(2)))~=8)
    error(['The value(s) in the second input ''DATES'' must be' ...
        ' numerical date(s) in YYYYMMDD format.'])
else
     DATEstart = DATES(1); DATEend = DATES(2);
end
if nargin < 3
    % Temporal resolutions, min/60 = hrs
    dt = 3; %dt = dt./60;
    dtskewn = 30; %dtskewn = dtskewn./60;
    weighting = true;
end
if nargin < 4
    % Temporal resolutions, min/60 = hrs
    dtskewn = 30; %dtskewn = dtskewn./60;
    if ~isnumeric(dt) || int16(dt)~=(dt)
        error(['The 3rd input must a numerical scalar or vector'...
            ' specifying the temporal resolution in full minutes.'])
    end
    weighting = true;
elseif nargin == 4
    if ~isnumeric(dt) || int16(dt)~=(dt)
        error(['The 3rd input must a numerical scalar or vector'...
            ' specifying the temporal resolution in full minutes.'])
    end
    if ~isnumeric(dtskewn) || int16(dtskewn)~=(dtskewn)
        error(['The 4th input must a numerical scalar or vector'...
            ' specifying the temporal resolution for skewness in full minutes.'])
    end
    weighting = true;
end
if nargin == 5
    if not(islogical(weighting))
        error('The 4th input has to be logical ''true'' or ''false''.')
    end
    if isempty(dt)
        dt = 3;
    end
    if isempty(dtskewn)
        dtskewn = 30;
    end
elseif nargin > 5
    error('Too many inputs.')
end

% Set time resolution for skewness
dtskewn = [num2str(dtskewn) 'min'];

for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):...
        datenum(num2str(DATEend),'yyyymmdd')

    % Check input and ouput files
    thedate = datestr(DATEi,'yyyymmdd');
    DATE = str2double(thedate);
    
    % Get default and site/unit/period specific parameters
    C = getconfig(site,DATE);
    
    % Get list of wstats files
    [dir_wstats_in, wstats_files] = getHALOfileList(site,DATE,'product','wstats');
    % If no files for today, skip the day
    if isempty(wstats_files)
        fprintf('\nNo ''wstats'' files found for ''%s'' at ''%s'', skipping...\n',thedate,site)
        continue
    end
        
     % Get list of tke files
    [dir_tke_in, tke_files] = getHALOfileList(site,DATE,'product','TKE');
    % If no files for today, skip the day
    if isempty(tke_files)
        fprintf('\nNo ''TKE'' files found for ''%s'' at ''%s'', skipping...\n',thedate,site)
        continue
    end
    
     % Get list of wind shear files
    [dir_shear_in, shear_files] = getHALOfileList(site,DATE,'product','windshear');
    % If no files for today, skip the day
    if isempty(shear_files)
        fprintf('\nNo ''windshear'' files found for ''%s'' at ''%s'', skipping...\n',thedate,site)
        continue
    end
    
    % Get & check output path, can the file be written
    [dir_BLc_out,~] = getHALOfileList(site,DATE,'product','ABLclassification');    
    status = checkHALOpath(site,DATE,'product','ABLclassification');
    if isempty(status)
        fprintf('Cannot write BL classification product for the site %s and date %s.',site,num2str(DATE));
        continue;
    end
        
    %% Set thresholds
    th.heatflux = 0;
    th.heatflux_hi = 10;
    th.epsilon = -5; % log10, is it turbulent?
    th.epsilon_hi = -4; % log10, is it convective?
    th.windshear = 0.03;
    th.vert_velo = 0;
    th.cloud = -5; % log10
    
    %% Load low level jets - TBD
    
    %% Load supplementary data - TBD

    for it = 1:length(dt)
        
        % temporal resolutions
        dt_i = [num2str(dt(it)) 'min'];

        %% Find the top of the aerosol layer       
        % load signal and beta
        wstats = load_nc_struct(fullfile([dir_wstats_in '/' wstats_files{1}]));
        signal_it = wstats.(['signal_mean_' dt_i]); 
        beta_it = real(log10(wstats.(['beta_mean_' dt_i])));
        
        % Filter out noise
        beta_it(10*real(log10(signal_it-1))<-20 | isnan(signal_it)) = nan;
        beta_top = beta_it; beta_top(beta_top > th.cloud) = nan;
        
        
        % search for the aerosol layer top
        itop_beta_it = ones(size(beta_it));
        aero_top_beta_it = nan(size(beta_it,1),1);
        aero_layer_mask_it = zeros(size(beta_it));
        for ip = 1:size(beta_it,1)
            if all(isnan(beta_it(ip,:)),2), continue; end
            for jp = 4:size(beta_it,2)-2
                itop_beta_it(ip,jp) = not((isnan(beta_top(ip,jp)) | beta_top(ip,jp) == 0) &...
                    (all(isnan(beta_top(ip,jp:jp+2))) | all(beta_top(ip,jp:jp+2)==0)));
            end
            itop_beta_it(ip,[1:3,find(itop_beta_it(ip,:) == 0,1,'last'):end]) = 0;
            tmp = itop_beta_it(ip,:); tmp(1:3) = 1;
            if not(isempty(find(tmp == 1,1,'last')))
                aero_top_beta_it(ip) = find(tmp == 1,1,'last')-1;
            end
            if not(isnan(aero_top_beta_it(ip)))
                aero_layer_mask_it(ip,1:aero_top_beta_it(ip)-2) = 1;
            end
        end
        
        % Clean up
        aero_layer_mask_it(:,1:3) = 0;
        aero_top_beta_it(aero_top_beta_it<4) = nan;
        aero_top_beta_it = round(smooth(aero_top_beta_it,5));

        % place into the output struct
        BLclass.(['aerosol_layer_top_' dt_i]) = aero_top_beta_it;
        BLclass.(['aerosol_layer_mask_' dt_i]) = aero_layer_mask_it;
        
        %% Load TKE
        TKE = load_nc_struct(fullfile([dir_tke_in '/' tke_files{1}]));

        %% Load wind shear
        windshear = load_nc_struct(fullfile([dir_shear_in '/' shear_files{1}]));

        %% Load skewness and re-grid into the temporal resolution of interest
        if weighting
        skewn_it = wstats.(['radial_velocity_weighted_skewness_' dtskewn]);
        skewn_error_it = wstats.(['radial_velocity_weighted_skewness_error_' dtskewn]);
        else
        skewn_it = wstats.(['radial_velocity_skewness_' dtskewn]);
        skewn_error_it = wstats.(['radial_velocity_skewness_error_' dtskewn]);
        end
        skewn_it(skewn_error_it>1) = 0;
        skewn_fill_tmp = inpaint_nans(skewn_it,4);
        [Xr,Yr] = meshgrid(wstats.height, wstats.(['time_' dt_i]));
        [Xo,Yo] = meshgrid(wstats.height,wstats.(['time_' dtskewn]));
        skewn_interp = interp2(Xo,Yo,skewn_fill_tmp,Xr,Yr,'linear');
        skewn_interp_filled = inpaint_nans(skewn_interp,4);
        skewn_interp_smooth = windowSlider(skewn_interp_filled,[10,5],@nanmedian);
        skewn_interp_smooth(isnan(beta_it)) = nan;

        %% Estimate sun rise and sun set times
        sun_rise_set = suncycle(C.latitude,C.longitude,datenum(thedate,'yyyymmdd'));

        %% Cloud mask
        cloudmask = zeros(size(beta_it));
        cloudbit = get_droplet_bit_matine(wstats.height,wstats.(['beta_mean_' dt_i]),0);
        cloudmask(beta_it>th.cloud) = 1;
        for i = 1:size(cloudmask,1)
            ibit = find(cloudbit(i),1,'first');
            imask = find(cloudmask(i),1,'first');
            if not(isempty(ibit) | isempty(imask)) && not(ibit == 0 | imask == 0) && imask < ibit
                cloudmask(i,1:ibit) = 0;
            end
        end
        
        % Dilate cloudmask
        shift_m1 = [cloudmask(:,2:end) zeros(size(cloudmask,1),1)];
        shift_p1 = [zeros(size(cloudmask,1),1) cloudmask(:,1:end-1)];
        cloudmask_dilated = cloudmask + shift_m1 + shift_p1;
        cloudmask_dilated(cloudmask_dilated<0) = 1;
        cloudmask_dilated(:,1:3) = 0;
%         cloudmask_dilated(BLclass.(['aerosol_layer_mask_' dt_i]) == 0) = 0;
        cloudmask = logical(cloudmask_dilated);
        
        %% Mask for precipitation
        velo = wstats.(['radial_velocity_mean_' dt_i]);
        velo(~aero_layer_mask_it|cloudmask) = nan;
        [velo_mean,a11,b22] = windowSlider(velo,[3,7],@prctile,[],95);
        precip_mask = velo_mean < -1 & b22./a11 > .8;
        precip_bit = any(precip_mask,2);
        
        %% Smooth TKE
        if weighting
        epsilon = TKE.(['epsilon_w_' dt_i]);
        epsilon_log10 = real(log10(epsilon));
        epsilon_log10(TKE.(['epsilon_error_' dt_i]) > 2 |...
            10*real(log10(signal_it-1))<-20) = -7;
        else
        epsilon = TKE.(['epsilon_' dt_i]);
        epsilon_log10 = real(log10(epsilon));
        epsilon_log10(TKE.(['epsilon_error_' dt_i]) > 2 |...
            10*real(log10(signal_it-1))<-20) = -7;
        
        end
        epsilon_log10(isnan(beta_it)) = nan;
        nan_profiles = all(isnan(epsilon_log10),2);
        eps_smooth_nonans = nan(size(epsilon_log10));
        eps_smooth_nonans(logical(aero_layer_mask_it)|cloudmask) = epsilon_log10(logical(aero_layer_mask_it)|cloudmask);
        eps_smooth_nonans_filled = inpaint_nans(eps_smooth_nonans,4);
        eps_smooth_nonans_filled(not(aero_layer_mask_it)) = nan;
        eps_smooth_nonans_filled_smoooth = windowSlider(eps_smooth_nonans_filled,[3,3],@nanmedian);
        eps_filled = nan(size(beta_it));
        eps_filled(~nan_profiles,:) = eps_smooth_nonans_filled_smoooth(~nan_profiles,:);
        eps_filled(isnan(beta_it)) = nan;

        %% Calculate turbulence coupling
        fubarfield.(['coupled_' dt_i]) = calculateHALOturbulenceCoupling(...
            eps_filled,skewn_interp_smooth,cloudmask,th,sun_rise_set,wstats.time_3min);
        
        %% Create bitfield
        bitfield.(['bits_' dt_i]) = createHALObitfield(th, wstats.(['time_' dt_i]),wstats.(['beta_mean_' dt_i]),...
            BLclass.(['aerosol_layer_top_' dt_i]),[],eps_filled,windshear.(['vector_wind_shear_' dt_i]),...
            sun_rise_set,fubarfield.(['coupled_' dt_i]),precip_bit);

        %% Generate boundary layer classification product
        [product, product_attribute, product_dimensions] = ...
            createHALOatmBLclassificationMasks(wstats.(['time_' dt_i]),wstats.height,...
            bitfield.(['bits_' dt_i]),fubarfield.(['coupled_' dt_i]),dt_i);
        
        fnames = fieldnames(product);
        for ifn = 1:length(fnames)
            if not(strcmp(fnames{ifn},'height'))
                data.(fnames{ifn}) = product.(fnames{ifn});
                att.(fnames{ifn}) = product_attribute.(fnames{ifn});
            end
        end
        fnames_d = fieldnames(product_dimensions);
        for ifn_d = 1:length(fnames_d)
            if not(strcmp(fnames_d{ifn_d},'height'))
                dim.(fnames_d{ifn_d}) = product_dimensions.(fnames_d{ifn_d});
            end
        end
    end
    
    
    % latitude, longitude, altitude, elevation
    if isfield(wstats,'latitude')
        data.latitude = wstats.latitude;
        data.longitude = wstats.longitude;
        data.altitude = wstats.altitude;
    elseif isfield(C,'latitude')
        data.latitude = C.latitude;
        data.longitude = C.longitude;
        data.altitude = C.altitude_in_meters;
    else
        data.latitude = 0;
        data.longitude = 0;
        data.altitude = 0;
    end

    
    % latitude
    att.latitude = create_attributes(...
        {},...
        'Latitude of lidar', ...
        'degrees_north');
    att.latitude.standard_name = 'latitude';
    % longitude
    att.longitude = create_attributes(...
        {},...
        'Longitude of lidar', ...
        'degrees_east');
    att.longitude.standard_name = 'longitude';
    % altitude
    att.altitude = create_attributes(...
        {},...
        'Height of instrument above mean sea level', ...
        'm');
    
    data.height = wstats.height;
    
    % Add height attribute
    att.height = create_attributes(...
        {'height'},...
        'Height above ground', ...
        'm',...
        [],...
        ['Range*sin(elevation), assumes that the instrument is at' ...
        ' ground level! If not, add the height of the instrument' ...
        ' from the ground to the height variable.']);
    
    % Height dimensions
    dim.height = length(wstats.height);
    
    % Create global attributs
    att.global.Conventions = 'CF-1.0';
    att.global.system = C.system;
    att.global.location = C.location;
    att.global.source = C.source;
    att.global.institution = C.institution;
    att.global.title = C.title;
    thedate = num2str(DATE);
    att.global.day   = str2double(thedate(7:8));
    att.global.month = str2double(thedate(5:6));
    att.global.year  = str2double(thedate(1:4));
    current_date = datestr(now,'yyyy-mm-dd HH:MM:SS');
    att.global.history = [current_date ' - Created by ' C.user];
    
    % Order fields
    data = orderfields(data);
    att  = orderfields(att);
    
    % Write into new netcdf
    write_nc_struct(fullfile([dir_BLc_out '/' thedate '_' site ...
        '_halo-doppler-lidar_BL-classification.nc']), ...
        dim, data, att)
end
end


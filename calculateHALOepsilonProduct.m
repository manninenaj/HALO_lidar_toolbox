function calculateHALOepsilonProduct(site,DATES,windproduct,typeof,weighting)
%calculateHALOepsilonProduct calculates the dissipation rate of 
%turbulent kinetic energy directly from vertical velocity variance with
%temporal resolution the Halo wstats product is given, and writes the 
%results into a daily *.nc file.
%
% Usage:
% calculateHALOverticalepsilonproduct(site,DATES,windproduct,typeof)
% calculateHALOverticalepsilonproduct(site,DATES,windproduct,typeof,weighting)
%
% Inputs:
% -site          String, site name, e.g. site = 'kuopio'
% -DATES         Scalar or vector, numeric, e.g. DATES = 20170401
%                or DATES = [20170401 20170431]
% -windproduct   String, name of the wind prodcut, 'windvad', 'winddbs'
% -typeof        String, with vad winds specify elevation angle (0-90),
%                e.g. '60', or with dbs winds specify number of beams, e.g. 
%                '3beams', '4beams'
% -weighting     logical 'true' or 'false'
%
% Created 2018-01-18, modified last 2019-09-17
% Antti Manninen
% Finnish Meteorological Institute
% antti.manninen@fmi.fi

if nargin < 4
  error('''site'', ''DATES'', ''windproduct'', and ''typeof'' are required inputs!')
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
if ~ischar(windproduct) || (~strcmp('windvad',windproduct) && ~strcmp('winddbs',windproduct))
  error('The 3rd input ''windproduct'' must be a string and can be either ''windvad'' or ''winddbs''.')
end
if strcmp(windproduct,'windvad')
    if ~ischar(typeof) || length(typeof) ~= 2 || (~isempty(str2num(typeof)) && str2num(typeof)<0 || str2num(typeof)>90)
        error('For vad winds, the 4th input must be a string and no longer than 2 characters specifying the elevation angle 0-90 degrees.')
    end
elseif strcmp(windproduct,'winddbs')
    if not(any(strcmp(typeof,{'3beams','4beams','5beams'})))
        error('For dbs winds, the 4th input must be a string and be in the form: ''3beams'', ''4beams'', or ''5beams''.')
    end
end
if nargin < 5
    weighting = false;
elseif nargin == 5
    if not(islogical(weighting))
        error('The 4th input has to be logical ''true'' or ''false''.')
    end
elseif nargin > 6
    error('Too many inputs.')
end

% Use datenum to accommodate leap years etc.
for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):...
            datenum(num2str(DATEend),'yyyymmdd')

    % Convert date into required formats
    thedate = datestr(DATEi,'yyyymmdd');
    thedate_yd = datestr(DATEi-1,'yyyymmdd');
    thedate_tw = datestr(DATEi+1,'yyyymmdd');
    DATE = str2double(thedate);
    DATE_yd = str2double(thedate_yd);
    DATE_tw = str2double(thedate_tw);


    % Get default and site/unit/period specific parameters
    C = getconfig(site,DATE);

    % Get list of files
    [dir_wstats_in, wstats_files] = getHALOfileList(site,DATE,'product','wstats');
    % If no files for today, skip the day
    if isempty(wstats_files)
        fprintf('\nNo ''wstats'' files found for ''%s'' at ''%s'', skipping...\n',thedate,site)
        continue
    end
    
    switch windproduct
        case 'windvad'
            typeof1 = ['ele' typeof];
        case 'winddbs'
            typeof1 = typeof;
    end
    [dir_wind_in, wind_files_tday] = getHALOfileList(site,DATE,'product',windproduct,typeof1);
    % If no files for today, skip the day
    if isempty(wind_files_tday)
        fprintf('\nNo ''%s'' (elev = %s) files found for ''%s'' at ''%s'', skipping...\n',windproduct,typeof,thedate,site)
        continue
    end
    [~, wind_files_yday] = getHALOfileList(site,DATE_yd,'product',windproduct,typeof1);
    [~, wind_files_tmrw] = getHALOfileList(site,DATE_tw,'product',windproduct,typeof1);


    % Get delta time from 'co' data (if 'wstats' exist, 'co' exists...)
    [dir_co_in, co_files] = getHALOfileList(site,DATE,'calibrated','stare','co');
    if isempty(co_files)
        fprintf('\nNo ''stare'' ''co'' files found for ''%s'' at ''%s'', skipping...\n',thedate,site)
        continue
    end
    data_co = load_nc_struct(fullfile([dir_co_in '/' co_files{1}]),{'time'});
    dt_raw = median(diff(data_co.time));

    % Get & check output path, can the file be written? 
    [dir_epsilon_out,~] = getHALOfileList(site,DATE,'product','epsilon');    
    status = checkHALOpath(site,DATE,'product','epsilon');
    if isempty(status)
        fprintf('Cannot write product epsilon for the site %s and date %s.',site,num2str(DATE));
        continue;
    end

    % Load, assume only one *.nc file per day
    wstats = load_nc_struct(fullfile([dir_wstats_in '/' wstats_files{1}]));
    wind_tday = load_nc_struct(fullfile([dir_wind_in '/' wind_files_tday{1}]));

    % Create common attribues, fields, and dimensions
    [data,att,dim] = createORcopyCommonAttsDims(wstats,'product',C);

    if strcmp(windproduct,'winddbs')
        if not(isfield(wind_tday,'height'))
            wind_tday.height = wind_tday.range(:).*sind(nanmedian(wind_tday.elevation(:)));
        end
        if not(isfield(wind_tday,'u'))
            % Calculate u and v components
            wind_tday.u = abs(wind_tday.wind_speed) .* sin(pi.*wind_tday.wind_direction/180);% u = u';
        end
        if not(isfield(wind_tday,'v'))
            wind_tday.v = abs(wind_tday.wind_speed) .* cos(pi.*wind_tday.wind_direction/180);% v = v';
        end
    end
    if ~isempty(wind_files_yday)
        % load if exists
        wind_yday = load_nc_struct(fullfile([dir_wind_in '/' wind_files_yday{1}]));
    else
        % empty otherwise
        wind_yday = [];
    end
    if ~isempty(wind_files_tmrw)
        % load if exists
        wind_tmrw = load_nc_struct(fullfile([dir_wind_in '/' wind_files_tmrw{1}]));
    else 
        % empty otherwise
        wind_tmrw = [];
    end
        
    % Get field names and extract the variables of interest
    fnames = fieldnames(wstats);
    ifs = strmatch('radial_velocity_simple_variance', fnames);
    fnames_var = fnames(ifs);
    if isempty(fnames_var)
        warning(['Simple variance not calculated for ' thedate ' at ' site '. Try recalculating wstats.'])
        continue
    end
    ifs_false = strmatch('radial_velocity_simple_variance_error',fnames_var);
    fnames_var(ifs_false) = [];

    % The naming of this variable has changed over the iterations. For existing data all versions are checked.
    name_options = cellstr(strvcat('radial_velocity_instrumental_error_mean',...
                                   'radial_velocity_instrumental_precision_mean',...
                                   'radial_velocity_instrumental_uncertainty_mean'));
    switch find(ismember(name_options,fnames))
      case 1 % error
          ifs_e_mean = strmatch('radial_velocity_instrumental_error_mean', fnames);
          ifs_e_var = strmatch('radial_velocity_instrumental_error_variance', fnames);
      case 2 % precision
          ifs_e_mean = strmatch('radial_velocity_instrumental_precision_mean', fnames);
          ifs_e_var = strmatch('radial_velocity_instrumental_precision_variance', fnames);       
      case 3 % uncertainty
          ifs_e_mean = strmatch('radial_velocity_instrumental_uncertainty_mean', fnames);
          ifs_e_var = strmatch('radial_velocity_instrumental_uncertainty_variance', fnames);               
    end
    
    fnames_mean_e = fnames(ifs_e_mean);
    fnames_var_e = fnames(ifs_e_var);
    
    ifs_nsamples = strmatch('nsamples', fnames);
    fnames_nsamples = fnames(ifs_nsamples);
    
    if weighting
        ifs_w = strmatch('radial_velocity_weighted_variance', fnames);
        fnames_wvar = fnames(ifs_w);
        ifs_wfalse = strmatch('radial_velocity_weighted_variance_error',fnames_wvar);
        fnames_wvar(ifs_wfalse) = [];

        % The naming of this variable has changed over the iterations. For existing data all versions are checked.
        wname_options = cellstr(strvcat('radial_velocity_instrumental_error_mean',...
                                        'radial_velocity_instrumental_precision_mean',...
                                        'radial_velocity_instrumental_uncertainty_mean'));
        switch find(ismember(wname_options,fnames_wvar))
          case 1
              ifs_e_mean = strmatch('radial_velocity_instrumental_error_weighted_mean', fnames_wvar);
              ifs_e_var = strmatch('radial_velocity_instrumental_error_weighted_variance', fnames_wvar);
          case 2
              ifs_e_mean = strmatch('radial_velocity_instrumental_precision_weighted_mean', fnames_wvar);
              ifs_e_var = strmatch('radial_velocity_instrumental_precision_weighted_variance', fnames_wvar);       
          case 3
              ifs_e_mean = strmatch('radial_velocity_instrumental_uncertainty_weighted_mean', fnames_wvar);
              ifs_e_var = strmatch('radial_velocity_instrumental_uncertainty_weighted_variance', fnames_wvar);               
        end

        fnames_wmean_e = fnames(ifs_e_wmean);    
        fnames_wvar_e = fnames(ifs_e_wvar);
    end
    
    ifs_t = strmatch('time', fnames);
    fnames_time = fnames(ifs_t);

    fprintf('\nGenerating the Halo epsilon product.\n')
    
    % Just to check how many time resolutions there are
    for ii = 1:length(fnames_var)
        dt_hrs = median(diff((wstats.(fnames_time{ii}))));
        tres = num2str(dt_hrs*60); % time reso in string       

        fprintf(['epsilon: estimating epsilon at ' tres ' min resolution...'])

        % true variance = variance - noise variance
        true_variance = wstats.(fnames_var{ii}) - ...
            wstats.(fnames_var_e{ii});
        
        if weighting
            % true weighted variance = weighted variance - noise weighted variance
            if isfield(wstats,'radial_velocity_weighted_mean_3min')
                true_wvariance = wstats.(fnames_wvar{ii})-wstats.(fnames_wvar_e{ii});
            else
                error('Weighted wstats have not been calculated for %s at %s.',thedate,site)
            end
        end

        % Filter by using number of samples
        nsamples_in_wstats = wstats.(fnames_nsamples{ii});
        
        % Re-grid the winds
        [Xr,Yr] = meshgrid(wstats.height, wstats.(fnames_time{ii}));
        if not(isempty(wind_tmrw) || isempty(wind_yday))
            [Xo,Yo] = meshgrid(wind_tday.height,...
                [wind_yday.time(:)-24;...
                wind_tday.time(:);...
                wind_tmrw.time(:)+24]);
            u_tmp = interp2(Xo,Yo,...
                [wind_yday.u;...
                wind_tday.u;...
                wind_tmrw.u],Xr,Yr);
            v_tmp = interp2(Xo,Yo,...
                [wind_yday.v;...
                wind_tday.v;...
                wind_tmrw.v],Xr,Yr);
            if isfield(wind_tday,'wind_speed_error')
            ws_e_tmp = interp2(Xo,Yo,...
                [wind_yday.wind_speed_error;...
                wind_tday.wind_speed_error;...
                wind_tmrw.wind_speed_error],Xr,Yr);
            else
                ws_e_tmp = repmat(1.5,length(wstats.(fnames_time{ii})),length(wstats.height)); % m s-1
            end
        else
            [Xo,Yo] = meshgrid(wind_tday.height, wind_tday.time);
            u_tmp = interp2(Xo,Yo,wind_tday.u,Xr,Yr);
            v_tmp = interp2(Xo,Yo,wind_tday.v,Xr,Yr);
            if isfield(wind_tday,'wind_speed_error')
                ws_e_tmp = interp2(Xo,Yo,wind_tday.wind_speed_error,Xr,Yr);
            else
                ws_e_tmp = repmat(1.5,length(wstats.(fnames_time{ii})),length(wstats.height)); % m s-1
            end
            
        end
        
        % Try infilling some gaps
        u_smooth = medianfilter(u_tmp);
        v_smooth = medianfilter(v_tmp);
        ws_e_smooth = medianfilter(ws_e_tmp);
        [~, count] = medianfilter(isfinite(u_tmp));
 
        % kernel is [3 3], max count is 9
        u_smooth(count < 3) = NaN;
        ind_nm = isnan(u_tmp) & isfinite(u_smooth);
        u = u_tmp; v = v_tmp; ws_e = ws_e_tmp;
        u(ind_nm) = u_smooth(ind_nm);
        v(ind_nm) = v_smooth(ind_nm);
        ws_e(ind_nm) = ws_e_smooth(ind_nm);
        
        % pass nans
        u(isnan(wstats.(fnames_var{ii}))) = nan;
        v(isnan(wstats.(fnames_var{ii}))) = nan;
        ws_e(isnan(wstats.(fnames_var{ii}))) = nan;
        
        % Re-calculate wind speed
        wind = sqrt(u .* u + v .* v);
        
        % Calculate beam diameter
        beamdiameter = ones(length(wstats.(fnames_time{ii})),1) * transpose(wstats.height) .* (C.beamwidth*pi/180);
        
        % Calculate length scales
        L_ii = ((3600.*dt_hrs.*wind) + 2 .* (ones(length(wstats.(fnames_time{ii})),1) * transpose(wstats.height)) .* sin(beamdiameter./2));
        L1_ii = ((3600.*dt_raw.*wind) + 2 .* (ones(length(wstats.(fnames_time{ii})),1) * transpose(wstats.height)) .* sin(beamdiameter./2));
        
        %     L1(isnan(L1)) = 0; L2(isnan(L2)) = 0;
        true_variance = abs(true_variance);
        true_variance(isnan(true_variance) | true_variance <= 0) = nan;
        if weighting
        true_wvariance(isnan(true_variance) | true_wvariance <= 0) = nan;
        end
        % Epsilon with true variance, instrumental noise variance removed
        epsilon_tmp = real((2/(3.*.55)).^1.5 .* true_variance.^(3/2) .* 2.*pi .* ((L_ii.^(2/3) - L1_ii.^(2/3)).^(-3/2)));
        epsilon_error_tmp = real(3.*sqrt(sqrt(4./C.num_samples_gate_stare_co.*(wstats.(fnames_mean_e{ii}).^2)./(wstats.(fnames_var{ii})))) + ws_e ./ L_ii);
        if weighting
        % Epsilon with weighted variance, instrumental noise variance taken into account
        epsilon_w_tmp = real((2/(3.*.55)).^1.5 .* true_wvariance.^(3/2) .* 2.*pi .* ((L_ii.^(2/3) - L1_ii.^(2/3)).^(-3/2)));
        epsilon_w_error_tmp = real(3.*sqrt(sqrt(4./C.num_samples_gate_stare_co.*(wstats.(fnames_wmean_e{ii}).^2)./(wstats.(fnames_wvar{ii})))) + ws_e ./ L_ii);
        end
        % final clean
        epsilon_tmp(~isfinite(epsilon_tmp)) = nan;
        epsilon_error_tmp(isnan(epsilon_tmp)) = nan;
        if weighting
        epsilon_w_tmp(~isfinite(epsilon_w_tmp)) = NaN;
        epsilon_w_error_tmp(isnan(epsilon_w_tmp)) = nan;
        end
            
        % Arrange into struct for writing
        data.(['time_' tres 'min']) = wstats.(fnames_time{ii});
        data.(['epsilon_' tres 'min']) = epsilon_tmp;
        data.(['epsilon_error_' tres 'min']) = epsilon_error_tmp;
        if weighting
        data.(['epsilon_w_' tres 'min']) = epsilon_w_tmp;
        data.(['epsilon_w_error_' tres 'min']) = epsilon_w_error_tmp;
        end
        data.(['L_' tres 'min']) = L_ii;
        data.(['L1_' tres 'min']) = L1_ii;
        
        % Create attributes
        % time
        att.(['time_' tres 'min']) = create_attributes(...
            {['time_' tres 'min']},...
            'Decimal hours UTC', ...
            'Hours UTC',...
            [],...
            ['Discrete time steps, in ' tres ' min temporal resolution.']);
        att.(['time_' tres 'min']).axis = 'T';
        % epsilon
        att.(['epsilon_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'].'range'},...
            'Dissipation rate of turbulent kinetic energy',...
            {'m2 s-3',''},...
            C.missing_value,...
            'Estimated from radial velocity variance unbiased by random noise and sample size, otherwise follows OConnor et al. (2010).',...
            {[-6 -1], 'logarithmic'});
        % epsilon error
        att.(['epsilon_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'].'range'},...
            'Fractional error in the dissipation rate of turbulent kinetic energy (epsilon)',...
            'unitless',...
            C.missing_value,...
            'Absolute error in epsilon divided by epsilon, OConnor et al. (2010).',...
            {[0 300], 'linear'});
        if weighting
        % epsilon_w
        att.(['epsilon_w_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'].'range'},...
            'Dissipation rate of turbulent kinetic energy',...
            {'m2 s-3',''},...
            C.missing_value,...
            'Estimated from weighted radial velocity variance unbiased by random noise and sample size, otherwise follows OConnor et al. (2010).',...
            {[-6 -1], 'logarithmic'});
        % epsilon_w error
        att.(['epsilon_w_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'].'range'},...
            'Fractional error in the dissipation rate of turbulent kinetic energy (epsilon) estimated from weighted radial velocity variance',...
            'unitless',...
            C.missing_value,...
            'Absolute error in epsilon divided by epsilon. Epsilon estimated from weighted radial velocity variance, OConnor et al. (2010).',...
            {[0 300], 'linear'});
        end
        % L1
        att.(['L_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'].'range'},...
            'Length scale of the largest eddies that pass completely through the lidar beam during the averaging window.',...
            {'m',''},...
            C.missing_value,...
            'Denoted as L in OConnor et al. (2010).',...
            {[0 4000], 'linear'});
        % L2
        att.(['L1_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'].'range'},...
            'Length scale of the scattering volume dimension per single sample in the averaging window.',...
            {'m',''},...
            C.missing_value,...
            'Outcome of Eq. (5) in OConnor et al. (2010).',...
            {[0 200], 'linear'});

        % Create dimensions
        dim.(['time_' tres 'min']) = length(data.(['time_' tres 'min']));
        fprintf('done.\n')    

    end   


    % Create global attributs
    att.global.Conventions = 'CF-1.0';
    att.global.system = C.system;
    att.global.location = C.location;
    att.global.source = C.source;
    att.global.institution = C.institution;
    att.global.title = C.title;
    att.global.day   = int16(str2double(thedate(7:8)));
    att.global.month = int16(str2double(thedate(5:6)));
    att.global.year  = int16(str2double(thedate(1:4)));
    current_date = datestr(now,'yyyy-mm-dd HH:MM:SS');
    att.global.history = [current_date ' - Created by ' C.user ];
    
    % Order fields
    data = orderfields(data);
    att  = orderfields(att);
    
    % Write into new netcdf
    write_nc_silent(fullfile([dir_epsilon_out '/' thedate ...
        '_' site '_halo-doppler-lidar_epsilon.nc']), dim, data, att)
end

end

    

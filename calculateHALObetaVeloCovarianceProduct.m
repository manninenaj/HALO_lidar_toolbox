function calculateHALObetaVeloCovarianceProduct(site,DATES,dt,drange)
%calculateHALObetaVeloCovarianceProduct calculates covariance between the
%attenuated backscatter and vertical velocity, which are read from the
%wstats product.
%
% Usage:
% calculateHALObetaVeloCovarianceProduct(site,DATES)
% calculateHALObetaVeloCovarianceProduct(site,DATES,dt)
% calculateHALObetaVeloCovarianceProduct(site,DATES,dt,drange)
%
% Inputs:
% -site        string, site name, e.g. site = 'kuopio'
% -DATES       scalar or vector, numeric, e.g. DATES = 20170401
%              or DATES = [20170401 20170431]
% -dt          scalar or vector, numeric, time window size in minutes,
%              e.g. dt = [10 30 60], or dt = 30. Default: dt = 30
% -drange      scalar or vector, numeric, range resolution in
%              number of range bins, e.g. drange = [1 6 10], or
%              drange = 6. Default: drange = 6
%
% Created 2018-10-28
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi


if nargin < 2
    error('''site'' and ''DATES'' are required inputs!')
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
    % temporal resolution, min / 3 min = window length
    dt = 90; dt = floor(dt./3);
    % range resolution, window heigth
    drange = 6;
elseif nargin == 3
    if ~isnumeric(dt) | int16(dt)~=dt | mod(dt,3)~=0
        error(sprintf(['The 3rd input must a numerical scalar or vector'...
            ' specifying the temporal resolution(s) in full minutes,'...
            ' \nand must be divisible by 3, which is the lowest ''wstats'''...
            'time resolution in minutes).']))
    else
        % temporal resolution, min / 3 min = window length
        dt = dt./3;
        % range resolution, window heigth
        drange = 6;
    end
end
if nargin < 4
    % range resolution, window heigth
    drange = 6;
end
if nargin == 4
    if ~isnumeric(dt) | int16(dt)~=dt | dt > 60
        error(sprintf(['The 3rd input must a numerical scalar or vector'...
            ' specifying the temporal resolution(s) in full minutes,'...
            ' \nand must be divisible by 3, which is the lowest ''wstats'''...
            'time resolution in minutes).']))
    else
        % temporal resolution, min / 3 min = window length
        dt = floor(dt./3);
    end
    if ~isnumeric(drange) | int16(drange)~=drange | drange > 10
        error(['The 4th input must a numerical scalar or vector'...
            ' specifying the range resolution(s) in no. of range bins,'...
            ' and be less than or equal to 10.'])
    end
end
if nargin > 4
    error('Too many inputs.')
end

% Use datenum to accommodate leap years etc.
for iDATE = datenum(num2str(DATEstart),'yyyymmdd'):...
        datenum(num2str(DATEend),'yyyymmdd')
    
    % Convert date into required formats
    thedate = datestr(iDATE,'yyyymmdd');
    DATE = str2double(thedate);
    
    % Get default and site/unit/period specific parameters
    C = getconfig(site,DATE);
    
    % Get paths to folders
    [dirto,files] = getHALOfileList(site,DATE,'product','wstats');
    [dir_to_folder_out,~] = getHALOfileList(site,DATE,'product','betavelocovariance');
    
    % Check path to write out
    status = checkHALOpath(site,DATE,'product','betavelocovariance');
    if isempty(status)
        fprintf('Can''t write %s - %s.',num2str(DATE),site);
        continue;
    end
    
    % Load, assume only one *.nc file per day
    if isempty(files), continue; end
    
    fprintf('\nGenerating the Halo beta-velo-covariance product.\n')
    
    % Load data
    wstats = load_nc_struct(fullfile([dirto files{1}]),{'time_3min',...
        'height','beta_mean_3min','signal_mean_3min','radial_velocity_mean_3min'});
    
    % Get beta and velocity
    beta = wstats.beta_mean_3min*1e6; % Mm-1 sr-1 m s-1
    velo = wstats.radial_velocity_mean_3min; % m
    
    % Clean up noise, 1.003 with 3 min averagin ok?
%     beta(wstats.signal_mean_3min<1.003) = nan; beta(:,1:3) = nan;
%     velo(wstats.signal_mean_3min<1.003) = nan; velo(:,1:3) = nan;
    
    % Set up the little-bag-of-bootstraps.
    lbob.subsample_size = .67;
    lbob.n_subsamples = 20;
    lbob.n_trials = 2; % Integer
    lbob.score_func = []; % set up later, separately for each calculation
    lbob.agg_func = @nanmean;
    lbob.x = [];
    lbob.y = [];
    lbob.scores = []; % Estimates generated on each iteration of the outer loop.
    lbob.mean = []; % Final value calculated in the outer loop
    lbob.lo_bound = []; % lower confidence bound
    lbob.hi_bound = []; % higher confidence bound
    lbob.std = []; % standard error
    
    tres = '3';
    
    for i = 1:length(dt)
        for j = 1:length(drange)
            
            % Generate window
            win = [dt(i) drange(j)];
            winlen = num2str(dt(i)*3); % in minutes
            winheight = num2str(abs(wstats.height(1)-wstats.height(drange(j)) ))  % in metres
            % Padding with nans
            beta_tmp = [nan(floor(win(1)/2),size(beta,2)); beta; nan(floor(win(1)/2),size(beta,2))];
            beta_padded = [nan(size(beta_tmp,1),floor(win(2)/2)), beta_tmp, nan(size(beta_tmp,1),floor(win(2)/2))];
            velo_tmp = [nan(floor(win(1)/2),size(velo,2)); velo; nan(floor(win(1)/2),size(velo,2))];
            velo_padded = [nan(size(velo_tmp,1),floor(win(2)/2)), velo_tmp, nan(size(velo_tmp,1),floor(win(2)/2))];
            
            % Construct indices of the running window
            a0 = nan(1,size(beta ,2)); a0(1) = 0; for ii = 2:size(beta,2), a0(ii) = a0(ii-1)+size(beta_padded,1); end
            rtmp1 = repmat(a0,1,size(beta,1));
            rtmp1s = sort(rtmp1);
            a1 = repmat(rtmp1s,win(2)*win(1),1);
            rtmp2 = repmat(1:size(beta,1),1,size(beta,2));
            a2 = repmat(rtmp2,win(2)*win(1),1);
            a3 = a1 + a2;
            rtmp = repmat(transpose(0:win(1)-1),1,size(beta,1)*size(beta,2));
            a4 = a3 + repmat(rtmp,win(2),1);
            a5 = repmat(sort(repmat(transpose((1:win(2))-1)*size(beta_padded,1),win(1),1)),1,size(beta,1)*size(beta,2));
            ind_running = a4 + a5;
            
            % Calculate number of total samples and number of non-nan samples
            nsamples = reshape(sum(~isnan(beta_padded(ind_running))),size(beta));
            
            fprintf(['\n  covariance of beta and Doppler velocity, ' ...
                    num2str(dt(i)) ' min resolution...'])
            X = beta_padded(ind_running); % assign X as beta
            Y = velo_padded(ind_running); % assign Y as vertical velocity
            Y_errors = zeros(size(ind_running)); % as dummy
            lbob.score_func = 'covariance'; % re-assign scoring function
            lbob = littleBagOfBootstraps(lbob,X,Y,Y_errors);
            fprintf('done.\n')
            
            lbob.best_estimate(lbob.best_estimate == 0) = nan;
            lbob.standard_error(lbob.standard_error == 0) = nan;
            lbob.confidence_interval_low(lbob.confidence_interval_low == 0) = nan;
            lbob.confidence_interval_high(lbob.confidence_interval_high == 0) = nan;
            
            % Reshape and clean up noise
            cond_clean = isnan(beta);
            cova_flux = reshape(lbob.best_estimate,size(beta));
            cova_flux(cond_clean) = nan;
            cova_flux_e = reshape(lbob.standard_error,size(beta));
            cova_flux_e(cond_clean) = nan;
            cova_flux_conf_int_lo = reshape(lbob.confidence_interval_low,size(beta));
            cova_flux_conf_int_lo(cond_clean) = nan;
            cova_flux_conf_int_hi = reshape(lbob.confidence_interval_high,size(beta));
            cova_flux_conf_int_hi(cond_clean) = nan;
            nsamples(cond_clean) = nan;
            
            %%-- Create output data struct --%%
            data.(['covariance_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = cova_flux;
            data.(['covariance_error_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = cova_flux_e;
            data.(['covariance_conf_int_lo_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = cova_flux_conf_int_lo;
            data.(['covariance_conf_int_hi_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = cova_flux_conf_int_hi;
            data.(['nsamples_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = nsamples;

%            data.(['covarflux']) = cova_flux;
%            data.(['covarflux_error']) = cova_flux_e;
%            data.(['covar_conf_int_lo']) = cova_flux_conf_int_lo;
%            data.(['covar_conf_int_hi']) = cova_flux_conf_int_hi;
%            data.(['nsamples']) = nsamples;

            
            %%-- Create attributes --%%
            % cov
           att.(['covariance_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = create_attributes(...
                {['time_' tres 'min'],'height'},...
                'Covariance of attenuated backscatter and vertical Doppler velocity',...
                {'Mm-1 sr-1 m s-1',''},...
                C.missing_value,...
                'Eq (1) in Engelmann et al., (2008).',...
                {[-0.1 0.1], 'linear'});
            % cov error
            att.(['covariance_error_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = create_attributes(...
                {['time_' tres 'min'],'height'},...
                'Standard error of covariance of attenuated backscatter and vertical Doppler velocity',...
                {'Mm-1 sr-1 m s-1',''},...
                C.missing_value,...
                'Eq (1) in Engelmann et al., (2008).',...
                {[-0.1 0.1], 'linear'});
            % cov conf int lo
            att.(['covariance_conf_int_lo_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = create_attributes(...
                {['time_' tres 'min'],'height'},...
                'Lower confidence interval of covariance of attenuated backscatter and vertical Doppler velocity',...
                {'Mm-1 sr-1 m s-1',''},...
                C.missing_value,...
                'Eq (1) in Engelmann et al., (2008).',...
                {[-0.1 0.1], 'linear'});
            % cov conf int hi
            att.(['covariance_conf_int_hi_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = create_attributes(...
                {['time_' tres 'min'],'height'},...
                'Higher confidence interval of covariance of attenuated backscatter and vertical Doppler velocity',...
                {'Mm-1 sr-1 m s-1',''},...
                C.missing_value,...
                'Eq (1) in Engelmann et al., (2008).',...
                {[-0.1 0.1], 'linear'});
            % nsamples
            att.(['nsamples_' tres 'min_window_' winlen 'min_by_' winheight 'm']) = create_attributes(...
                {['time_' tres 'min'],'height'},...
                'Number of non-nan samples used in the calculation of covariance',...
                {'Mm-1 sr-1 m s-1',''},...
                C.missing_value,...
                [],...
                {[0 1000], 'linear'});

        end
    end
    
    
    % time
    data.(['time_' tres 'min']) = wstats.(['time_' tres 'min']);
    att.(['time_' tres 'min']) = create_attributes(...
        {['time_' tres 'min']},...
        'Decimal hours UTC', ...
        'Hours UTC',...
        [],...
        ['Running window, time of center point, in ' tres ' min temporal resolution.']);
    att.(['time_' tres 'min']).axis = 'T';
    
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
    
    % Height
    data.height = wstats.height;
    att.height = create_attributes(...
        {'height'},...
        'Height above ground', ...
        'm',...
        [],...
        ['Range*sin(elevation), assumes lidar is at ground level.']);
    
    % Add height dim
    dim.height = length(data.height);
    dim.(['time_' tres 'min']) = length(data.(['time_' tres 'min']));

    % Create global attributs
    att.global.Conventions = 'CF-1.0';
    att.global.system = C.system;
    att.global.location = C.location;
    att.global.source = C.source;
    att.global.institution = C.institution;
    att.global.title = C.title;
    att.global.day   = str2double(thedate(7:8));
    att.global.month = str2double(thedate(5:6));
    att.global.year  = str2double(thedate(1:4));
    current_date = datestr(now,'yyyy-mm-dd HH:MM:SS');
    att.global.history = [current_date ' - Created by ' C.user ];
    
    % Order fields
    data = orderfields(data);
    att  = orderfields(att);

data
att
dim
    
    % Write into new netcdf
%    write_nc_silent(fullfile([dir_to_folder_out '/' thedate ...
%        '_' site '_halo-doppler-lidar_covariance-beta-velo.nc']), dim, data, att)
 

fname = fullfile([dir_to_folder_out '/' thedate ...
        '_' site '_halo-doppler-lidar_covariance-beta-velo.nc'])

write_nc_struct(fullfile([dir_to_folder_out '/' thedate ...
        '_' site '_halo-doppler-lidar_covariance-beta-velo.nc']), dim, data, att)
write_nc_struct(['/home/cloudnet/' thedate ...
        '_' site '_halo-doppler-lidar_covariance-beta-velo.nc'], dim, data, att)

fname1 = ['/home/cloudnet/' thedate ...
        '_' site '_halo-doppler-lidar_covariance-beta-velo.nc']


[p a d]=load_nc_struct(fname1)
write_nc_struct(fname,d, p, a)                                                     


    
end
% % %     % Calculate covariance
% % %     [Cova_flux,~,nsamples] = ...
% % %         windowSlider2D(beta*1e6,[30,6],@covariance,'nans',[],velo);


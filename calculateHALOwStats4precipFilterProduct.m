function calculateHALOwStats4precipFilterProduct(site,DATES,dt,timeStep)
%calculateHALOwStatsProduct reads calibrated *co.nc files and
%calculates vertical radial velocity statistics: mean, std.dev.,
%variance, skewness, kurtosis with a discrete time steps. The results
%are written in a daily *.nc files. Also, the theoretical measurement
%errors, i.e. noise std. dev. and variance are calculated.
%
% Usage:
% calculateHALOwStatsProduct(site,DATES)
% calculateHALOwStatsProduct(site,DATES,dt)
% calculateHALOwStatsProduct(site,DATES,dt,weighting)
%
% Inputs:
% -site        string, site name, e.g. site = 'kuopio'
% -DATES       scalar or vector, numeric, e.g. DATES = 20170401
%              or DATES = [20170401 20170431]
% -dt          scalar or vector, numeric, temporal resolution in minutes,
%              e.g. dt = [1 10], or dt = 30. Default: dt = [3 10 30]
% -timeStep    specifies in minutes how long chunks are loaded and
%              processed at a time. Default: 240
% -weighting   logical 'true' or 'false'. Default: 'true'
%
% Created 2018-01-18
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
    % Temporal resolutions, min/60 = hrs
    dt = [3 30 60]; dt = dt./60;
    timeStep = 120;
elseif nargin == 3
    if ~isnumeric(dt) | int16(dt)~=dt | dt > 60
        error(['The 3rd input must a numerical scalar or vector'...
            ' specifying the temporal resolution in full minutes,'...
            ' and less than or equal to 60.'])
    else
        % Temporal resolutions, min/60 = hrs
        dt = dt./60;
    end
end
if nargin < 4
    timeStep = 120;
end
if nargin == 4
    if ~isnumeric(dt) | int16(dt)~=dt | dt > 60
        error(['The 3rd input must a numerical scalar or vector'...
            ' specifying the temporal resolution in full minutes,'...
            ' and less than or equal to 60.'])
    else
        % Temporal resolutions, min/60 = hrs
        dt = dt./60;
    end
    if ~isnumeric(timeStep) | timeStep < 60 | timeStep > 1440 | ~isscalar(timeStep) | ...
            int16(timeStep)~=timeStep | ~all(int16(timeStep./dt)==timeStep./dt)
        error(sprintf(['The 4th input ''timeSptep'' has to be numerical integer,'...
            ' larger than 1 (min), less than 1440 (min),\nand'...
            ' divisible by the 3rd input ''dt''. So, function works only'...
            ' with multiples of 60 (min): 60, 120, 240.']))
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
    
    % Get list of files
    [dir_to_folder_in,halo_files] = getHALOfileList(site,DATE,'calibrated','stare','co');
    [dir_to_folder_out,~] = getHALOfileList(site,DATE,'product','wstats4precipfilter');
    
    % Check path to write out
    status = checkHALOpath(site,DATE,'product','wstats4precipfilter');
    if isempty(status)
        fprintf('Can''t write %s - %s.',num2str(DATE),site);
        continue;
    end
    
    % Load, assume only one *_co.nc file per day
    if isempty(halo_files), continue; end
    
    fprintf('\nGenerating the Halo wstats-4-precip-filter product.\n')
    
    % Define groups, which are used to process the data in chunks
    ncid = netcdf.open([dir_to_folder_in halo_files{1}],'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid,'time');
    time_info = double(netcdf.getVar(ncid,varid));
    varid = netcdf.inqVarID(ncid,'range');
    range_len = length(double(netcdf.getVar(ncid,varid)));
    netcdf.close(ncid)
    
    ref_time_info = transpose((1:1:86400)/3600);
    ref_time_block_indices = ones(timeStep*60,86400/(timeStep*60)) + ...
        repmat(0:(86400/(timeStep*60))-1,timeStep*60,1);
    ref_time_block_indices = ref_time_block_indices(:);
    [~,ib] = look4nearest(time_info,ref_time_info);
    ref_time_info_block = ref_time_block_indices(ib);
    
    att = []; dim = [];
    for idt1 = 1:length(dt)
        tres = num2str(dt(idt1)*60); % time reso in string
        bn = 'radial_velocity'; % base for the name
        
        att = addAttributes(C,tres,bn,att);
        
        % time
        data.(['time_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        dim.(['time_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        
        % Initialize cell arrays for collection of the chunks
        data.([bn '_mean_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.([bn '_mean_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_mean_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_variance_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_mean_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_variance_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_mean_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_variance_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_mean_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_variance_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);

        data.([bn '_wmean_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.([bn '_wmean_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_wmean_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_wvariance_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_wmean_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_wvariance_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_wmean_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_wvariance_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_wmean_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['beta_wvariance_error_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);

        
        data.(['signal_instrumental_precision_mean_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['signal_instrumental_precision_variance_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
        data.(['nsamples_' tres 'min']) = cell(numel(unique(ref_time_block_indices)),1);
    end
    
    % Iterate over unique time steps
    for ichunk = 1:numel(unique(ref_time_block_indices))
        
        % Find start and end of the time window
        ibegin = find(ref_time_info_block == ichunk,1,'first');
        iend = find(ref_time_info_block == ichunk,1,'last');
        if isempty(ibegin) || isempty(iend), continue; end
        
        fprintf('\nwstats: loading %s and calculating block %s/%s ...\n', ...
            datestr(iDATE,'yyyymmdd'),num2str(ichunk), ...
            num2str(numel(unique(ref_time_block_indices))))
        
        % Extract data
        istarts = [0,ibegin-1];
        icounts = [range_len,iend-ibegin];
        ncid = netcdf.open([dir_to_folder_in halo_files{1}],'NC_NOWRITE');
        varid = netcdf.inqVarID(ncid,'v_raw');
        tmp.v_raw = transpose(double(netcdf.getVar(ncid,varid,istarts,icounts)));
        varid = netcdf.inqVarID(ncid,'v_error');
        tmp.v_error = transpose(double(netcdf.getVar(ncid,varid,istarts,icounts)));
        varid = netcdf.inqVarID(ncid,'signal');
        tmp.signal = transpose(double(netcdf.getVar(ncid,varid,istarts,icounts)));
        varid = netcdf.inqVarID(ncid,'beta_raw');
        tmp.beta_raw = transpose(double(netcdf.getVar(ncid,varid,istarts,icounts)));
        varid = netcdf.inqVarID(ncid,'beta_error');
        tmp.beta_error = transpose(double(netcdf.getVar(ncid,varid,istarts,icounts)));
        varid = netcdf.inqVarID(ncid,'time');
        tmp.time = double(netcdf.getVar(ncid,varid,istarts(2),icounts(2))); tmp.time = tmp.time(:);
        varid = netcdf.inqVarID(ncid,'range');
        tmp.range = double(netcdf.getVar(ncid,varid,istarts(1),icounts(1))); tmp.range = tmp.range(:);
        netcdf.close(ncid)
        
        % Check and equalize nans
        cond_nan = isnan(tmp.signal) | isnan(tmp.beta_raw) | isnan(tmp.v_raw) | ...
            tmp.signal == C.missing_value | tmp.beta_raw == C.missing_value | ...
            tmp.v_raw == C.missing_value;
        tmp.v_raw(cond_nan | tmp.v_raw == C.missing_value) = nan;
        tmp.v_error(cond_nan | tmp.v_error == C.missing_value) = nan;
        tmp.signal(cond_nan | tmp.signal == C.missing_value) = nan;
        tmp.beta_raw(cond_nan | tmp.beta_raw == C.missing_value) = nan;
        tmp.beta_error(cond_nan | tmp.beta_error == C.missing_value) = nan;
            
        % Grid into 1 sec resolution
        ref_time_sec = ref_time_info(ref_time_block_indices==ichunk);
        [~,ib] = look4nearest(tmp.time,ref_time_sec);
        ref_v_raw = nan(length(ref_time_sec),length(tmp.range));
        ref_v_error = nan(length(ref_time_sec),length(tmp.range));
        ref_signal = nan(length(ref_time_sec),length(tmp.range));
        ref_beta = nan(length(ref_time_sec),length(tmp.range));
        ref_beta_error = nan(length(ref_time_sec),length(tmp.range));
        ref_v_raw(ib,:) = tmp.v_raw;
        ref_v_error(ib,:) = tmp.v_error;
        ref_signal(ib,:) = tmp.signal;
        ref_beta(ib,:) = tmp.beta_raw;
        ref_beta_error(ib,:) = tmp.beta_error;
                
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
        
        % Iterate over the selected temporal resolutions
        for idt = 1:length(dt)
            
            tres = num2str(dt(idt)*60); % time reso in string
            bn = 'radial_velocity'; % base for the name
            
            % Create time steps
            atime = transpose(dt(idt)/2:dt(idt):24 - dt(idt)/2);
            
            % Create chunkt time refence vector.
            chunktime = sort(repmat(1:length(1:timeStep/60:24),1,24/dt(idt)/length(1:timeStep/60:24)));
            
            fprintf('\nwstats: %s min resolution.', num2str(dt(idt).*60))
            % Vectorize indices
            iv = reshape(1:numel(ref_v_raw),(dt(idt)*3600),size(ref_v_raw,1)/(dt(idt)*3600)*size(ref_v_raw,2));
            data.(['nsamples_' tres 'min']){ichunk} = reshape(sum(~isnan(ref_v_raw(iv))),size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            
            %%--- UNWEIGHTED STATISTICS ---%%
            fprintf('\n  radial velocity unweighted mean ...')
            X = []; % assign X as empty
            Y = ref_v_raw(iv); % assign Y as vertical velocity
            Y_errors = ref_v_error(iv); % assign Y_errors 
            lbob.score_func = @weightedMean; % re-assign scoring function
            lbob_velo_mean = littleBagOfBootstraps(lbob,X,Y,Y_errors,'unweighted',false);
            fprintf('done.')
            
            fprintf('\n  signal unweighted mean...')
            X = []; % assign X as empty
            Y = ref_signal(iv); % re-assign Y as signal
            Y_errors = ref_signal(iv).* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedMean;
            lbob_signal_mean = littleBagOfBootstraps(lbob,X,Y,Y_errors,'unweighted',false);
            fprintf('done.')
            
            fprintf('\n  signal unweighted variance...')
            X = []; % assign X as empty
            Y = ref_signal(iv); % re-assign Y as signal
            Y_errors = ref_signal(iv).* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedVariance;
            lbob_signal_var = littleBagOfBootstraps(lbob,X,Y,Y_errors,'unweighted',false);
            fprintf('done.')
            
            fprintf('\n  beta unweighted mean...')
            X = []; % assign X as empty
            Y = ref_beta(iv); % re-assign Y as beta
            Y_errors = ref_beta(iv).* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedMean;
            lbob_beta_mean = littleBagOfBootstraps(lbob,X,Y,Y_errors,'unweighted',false);
            fprintf('done.')
            
            fprintf('\n  beta unweighted variance...')
            X = []; % assign X as empty
            Y = ref_beta(iv); % re-assign Y as beta
            Y_errors = ref_beta(iv) .* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedVariance;
            lbob_beta_var = littleBagOfBootstraps(lbob,X,Y,Y_errors,'unweighted',false);
            fprintf('done.')
            
             %%--- WEIGHTED STATISTICS ---%%
            fprintf('\n  radial velocity weighted mean ...')
            X = []; % assign X as empty
            Y = ref_v_raw(iv); % assign Y as vertical velocity
            Y_errors = ref_v_error(iv); % assign Y_errors 
            lbob.score_func = @weightedMean; % re-assign scoring function
            lbob_velo_wmean = littleBagOfBootstraps(lbob,X,Y,Y_errors,'weighted',false);
            fprintf('done.')
            
            fprintf('\n  signal weighted mean...')
            X = []; % assign X as empty
            Y = ref_signal(iv); % re-assign Y as signal
            Y_errors = ref_signal(iv).* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedMean;
            lbob_signal_wmean = littleBagOfBootstraps(lbob,X,Y,Y_errors,'weighted',false);
            fprintf('done.')
            
            fprintf('\n  signal weighted variance...')
            X = []; % assign X as empty
            Y = ref_signal(iv); % re-assign Y as signal
            Y_errors = ref_signal(iv).* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedVariance;
            lbob_signal_wvar = littleBagOfBootstraps(lbob,X,Y,Y_errors,'weighted',false);
            fprintf('done.')
            
            fprintf('\n  beta weighted mean...')
            X = []; % assign X as empty
            Y = ref_beta(iv); % re-assign Y as beta
            Y_errors = ref_beta(iv).* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedMean;
            lbob_beta_wmean = littleBagOfBootstraps(lbob,X,Y,Y_errors,'weighted',false);
            fprintf('done.')
            
            fprintf('\n  beta weighted variance...')
            X = []; % assign X as empty
            Y = ref_beta(iv); % re-assign Y as beta
            Y_errors = ref_beta(iv) .* ref_beta_error(iv); % beta_error is fractional error
            lbob.score_func = @weightedVariance;
            lbob_beta_wvar = littleBagOfBootstraps(lbob,X,Y,Y_errors,'weighted',false);
            fprintf('done.')
            
            
             %%--- THE REST ---%%
            fprintf('\n  signal instrumental precision unweighted mean...')
            X = []; % assign X as empty
            Y = ref_beta_error(iv); % re-assign Y as beta
            lbob_signal_instrumental_precision_mean = nanmean(Y);
            fprintf('done.')
            
            fprintf('\n  signal instrumental precision unweighted variance...')
            X = []; % assign X as empty
            Y = ref_beta_error(iv); % re-assign Y as beta
            lbob_signal_instrumental_precision_var = nanvar(Y);
            fprintf('done.\n')
            
            % wstats
            lbob_velo_mean.best_estimate(lbob_velo_mean.best_estimate == 0) = nan;
            % wstats standard errors
            lbob_velo_mean.standard_error(lbob_velo_mean.standard_error == 0) = nan;
            % signal
            lbob_signal_mean.best_estimate(lbob_signal_mean.best_estimate == 0) = nan;
            lbob_signal_var.best_estimate(lbob_signal_var.best_estimate == 0) = nan;
            lbob_signal_mean.standard_error(lbob_signal_mean.standard_error == 0) = nan;
            lbob_signal_var.standard_error(lbob_signal_var.standard_error == 0) = nan;
            % beta
            lbob_beta_mean.best_estimate(lbob_beta_mean.best_estimate == 0) = nan;
            lbob_beta_var.best_estimate(lbob_beta_var.best_estimate == 0) = nan;
            lbob_beta_mean.standard_error(lbob_beta_mean.standard_error == 0) = nan;
            lbob_beta_var.standard_error(lbob_beta_var.standard_error == 0) = nan;

            % w wstats
            lbob_velo_wmean.best_estimate(lbob_velo_wmean.best_estimate == 0) = nan;
            % wstats standard errors
            lbob_velo_wmean.standard_error(lbob_velo_wmean.standard_error == 0) = nan;
            % w signal
            lbob_signal_wmean.best_estimate(lbob_signal_wmean.best_estimate == 0) = nan;
            lbob_signal_wvar.best_estimate(lbob_signal_wvar.best_estimate == 0) = nan;
            lbob_signal_wmean.standard_error(lbob_signal_wmean.standard_error == 0) = nan;
            lbob_signal_wvar.standard_error(lbob_signal_wvar.standard_error == 0) = nan;
            % w beta
            lbob_beta_wmean.best_estimate(lbob_beta_wmean.best_estimate == 0) = nan;
            lbob_beta_wvar.best_estimate(lbob_beta_wvar.best_estimate == 0) = nan;
            lbob_beta_wmean.standard_error(lbob_beta_wmean.standard_error == 0) = nan;
            lbob_beta_wvar.standard_error(lbob_beta_wvar.standard_error == 0) = nan;

            
            % beta noise
            lbob_signal_instrumental_precision_mean(lbob_signal_instrumental_precision_mean == 0) = nan;
            lbob_signal_instrumental_precision_var(lbob_signal_instrumental_precision_var == 0) = nan;
            
            % Collect
            data.(['time_' tres 'min']){ichunk} = atime(chunktime == ichunk);
            dim.(['time_' tres 'min']){ichunk} = length(atime(chunktime == ichunk));
            
            data.([bn '_mean_' tres 'min']){ichunk,1} = reshape(lbob_velo_mean.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.([bn '_mean_error_' tres 'min']){ichunk} = reshape(lbob_velo_mean.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            
            data.(['signal_mean_' tres 'min']){ichunk} = reshape(lbob_signal_mean.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['signal_variance_' tres 'min']){ichunk} = reshape(lbob_signal_var.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_mean_' tres 'min']){ichunk} = reshape(lbob_beta_mean.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_variance_' tres 'min']){ichunk} = reshape(lbob_beta_var.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            
            data.(['signal_mean_error_' tres 'min']){ichunk} = reshape(lbob_signal_mean.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['signal_variance_error_' tres 'min']){ichunk} = reshape(lbob_signal_var.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_mean_error_' tres 'min']){ichunk} = reshape(lbob_beta_mean.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_variance_error_' tres 'min']){ichunk} = reshape(lbob_beta_var.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));

            data.([bn '_wmean_' tres 'min']){ichunk,1} = reshape(lbob_velo_wmean.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.([bn '_wmean_error_' tres 'min']){ichunk} = reshape(lbob_velo_wmean.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            
            data.(['signal_wmean_' tres 'min']){ichunk} = reshape(lbob_signal_wmean.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['signal_wvariance_' tres 'min']){ichunk} = reshape(lbob_signal_wvar.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_wmean_' tres 'min']){ichunk} = reshape(lbob_beta_wmean.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_wvariance_' tres 'min']){ichunk} = reshape(lbob_beta_wvar.best_estimate,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            
            data.(['signal_wmean_error_' tres 'min']){ichunk} = reshape(lbob_signal_wmean.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['signal_wvariance_error_' tres 'min']){ichunk} = reshape(lbob_signal_wvar.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_wmean_error_' tres 'min']){ichunk} = reshape(lbob_beta_wmean.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['beta_wvariance_error_' tres 'min']){ichunk} = reshape(lbob_beta_wvar.standard_error,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));     
            
            data.(['signal_instrumental_precision_mean_' tres 'min']){ichunk} = reshape(lbob_signal_instrumental_precision_mean,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            data.(['signal_instrumental_precision_variance_' tres 'min']){ichunk} = reshape(lbob_signal_instrumental_precision_var,size(ref_v_raw,1)/(dt(idt)*3600),size(ref_v_raw,2));
            
        end
    end
    
    % Convert from cell to mat
    fnames = fieldnames(data);
    for ifn = 1:length(fnames)
        if iscell(data.(fnames{ifn}))
            data.(fnames{ifn}) = cell2mat(data.(fnames{ifn}));
            if strcmp(fnames{ifn},['time_' num2str(dt*60) 'min'])
                data.(fnames{ifn}) = data.(fnames{ifn})(:);
            end
        end
    end
    fnames = fieldnames(dim);
    for ifn = 1:length(fnames)
        if iscell(dim.(fnames{ifn}))
            dim.(fnames{ifn}) = sum(cell2mat(dim.(fnames{ifn})));
        end
    end
    %%-- Add variables --%%
    data.height = tmp.range;
    
    % latitude, longitude, altitude, elevation
    if isfield(tmp,'latitude')
        data.latitude = tmp.latitude;
        data.longitude = tmp.longitude;
        data.altitude = tmp.altitude;
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
    
    % height
    att.height = create_attributes(...
        {'height'},...
        'Height above ground', ...
        'm',...
        [],...
        ['Range*sin(elevation), assumes that the instrument is at ground level! If'...
        ' not, add the height of the instrument from the ground to the height variable.']);
    
    % Add height dim
    dim.height = length(data.height);
    
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
    
    % Order fivlds
    data = orderfields(data);
    att  = orderfields(att);
    
    % Write into new netcdf
    write_nc_struct(fullfile([dir_to_folder_out '/' thedate ...
        '_' site '_halo-doppler-lidar_wstats-4-precip-filter.nc']), dim, data, att)
end

    function [att] = addAttributes(C,tres,bn,att)
        
        %%-- ATTRIBUTES --%%
        % time
        att.(['time_' tres 'min']) = create_attributes(...
            {['time_' tres 'min']},...
            'Decimal hours UTC', ...
            'Hours UTC',...
            [],...
            ['Discrete time steps, in ' tres ' min temporal resolution.']);
        att.(['time_' tres 'min']).axis = 'T';
        
        % Variables
        % nsamples
        att.(['nsamples_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Number of samples (' tres ' min)'],...
            {'count','count'},...
            C.missing_value,...
            'Number of non-nans');  
        
        %%%---unweighted---%%
        % velo mean
        att.([bn '_mean_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Radial velocity mean (' tres ' min)'],...
            {'m s-1','m s<sup>-1</sup>'},...
            C.missing_value,...
            'Unbiased by random noise and sample size.',...
            {[-3 3], 'linear'});
        
        % Standard errors
        % velo mean error
        att.([bn '_mean_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in radial velocity mean (' tres ' min)'],...
            {'m s-1','m s<sup>-1</sup>'},...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});

        % Other variables
        % signal mean
        att.(['signal_mean_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Signal (raw) mean', ...
            'arbitrary units', ...
            C.missing_value, ...
            'Unbiased by random noise and sample size.');
        % signal variance
        att.(['signal_variance_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Signal (raw) variance', ...
            'arbitrary units', ...
            C.missing_value, ...
            'Unbiased by random noise and sample size.');
        % beta mean
        att.(['beta_mean_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Raw attenuated beta mean',...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            ['Unbiased by random noise and sample size.'],...
            {[1e-7 1e-4],'logarithmic'});
        % beta variance
        att.(['beta_variance_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Raw attenuated beta variance',...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            ['Unbiased by random noise and sample size.'],...
            {[1e-7 1e-4],'logarithmic'});
        
        % signal mean error
        att.(['signal_mean_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in mean of signal (' tres ' min)'],...
            'arbitrary units',...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        % signal var error
        att.(['signal_variance_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in variance of signal (' tres ' min)'],...
            'arbitrary units',...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        % beta mean error
        att.(['beta_mean_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in mean of beta (' tres ' min)'],...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        % beta var error
        att.(['beta_variance_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in variance of beta (' tres ' min)'],...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        
        
        %%%---weighted---%%
        % velo wmean
        att.([bn '_wmean_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Radial velocity mean (' tres ' min)'],...
            {'m s-1','m s<sup>-1</sup>'},...
            C.missing_value,...
            'Unbiased by random noise and sample size.',...
            {[-3 3], 'linear'});
        
        % Standard errors
        % velo wmean error
        att.([bn '_wmean_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in radial velocity mean (' tres ' min)'],...
            {'m s-1','m s<sup>-1</sup>'},...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});

        % Other variables
        % signal wmean
        att.(['signal_wmean_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Signal (raw) mean', ...
            'arbitrary units', ...
            C.missing_value, ...
            'Unbiased by random noise and sample size.');
        % signal wvariance
        att.(['signal_wvariance_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Signal (raw) variance', ...
            'arbitrary units', ...
            C.missing_value, ...
            'Unbiased by random noise and sample size.');
        % beta wmean
        att.(['beta_wmean_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Raw attenuated beta mean',...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            ['Unbiased by random noise and sample size.'],...
            {[1e-7 1e-4],'logarithmic'});
        % beta wvariance
        att.(['beta_wvariance_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            'Raw attenuated beta variance',...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            ['Unbiased by random noise and sample size.'],...
            {[1e-7 1e-4],'logarithmic'});
        
        % signal wmean error
        att.(['signal_wmean_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in mean of signal (' tres ' min)'],...
            'arbitrary units',...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        % signal wvar error
        att.(['signal_wvariance_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in variance of signal (' tres ' min)'],...
            'arbitrary units',...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        % betaw mean error
        att.(['beta_wmean_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in mean of beta (' tres ' min)'],...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        % beta wvar error
        att.(['beta_wvariance_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Standard error in variance of beta (' tres ' min)'],...
            {'sr-1 m-1','sr<sup>-1</sup> m<sup>-1</sup>'},...
            C.missing_value,...
            '',...
            {[0.01 10], 'logarithmic'});
        
        
        
        
        
        % signal instrumental precision mean
        att.(['signal_instrumental_precision_mean_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Mean of fractional instrumental precision in beta (' tres ' min)'],...
            {'m s-1','m s<sup>-1</sup>'},...
            C.missing_value,...
            'Mean of fractional error in beta. No correction for bias due to noise and sample size) correction.',...
            {[0 1], 'linear'});
        % signal instrumental precision variance
        att.(['signal_instrumental_precision_variance_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'height'},...
            ['Variance of fractional instrumental precision in beta (' tres ' min)'],...
            {'m s-1','m s<sup>-1</sup>'},...
            C.missing_value,...
            'Variance of fractional error in beta. No correction for bias due to noise and sample size) correction.',...
            {[0 1], 'log'});
        
    end
end



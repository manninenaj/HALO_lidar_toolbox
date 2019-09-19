function calculateHALOcloudProduct(site,DATES)
%calculateHALOcloudProduct calculates cloud base height, cloud base
%velocity, and provides cloud-precipitation mask in temporal resolution(s)
%based on what is/are available in the wstats product from the site and
%date in question.
%
% Usage:
% calculateHALOcloudProduct(site,DATES)
%
% Inputs:
% -site        string, site name, e.g. site = 'kuopio'
% -DATES       scalar or vector, numeric, e.g. DATES = 20170401 
%              or DATES = [20170401 20170431]
%
% Created 2018-01-18
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

% Check inputs
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

% Use datenum to accommodate leap years etc.
for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):...
        datenum(num2str(DATEend),'yyyymmdd')

    % Convert date into required formats
    thedate = datestr(DATEi,'yyyymmdd');
    DATE = str2double(thedate);
    
    % Get default and site/unit/period specific parameters
    C = getconfig(site,DATE);
       
    % Get list of files
        [dir_to_folder_in,halo_files] = getHALOfileList(site,DATE,'product','wstats');

        % Check path to write out
    status = checkHALOpath(site,DATE,'product','cloud');
    if isempty(status)
       fprintf('Can''t write %s - %s.',num2str(DATE),site);
       continue;
    end
    [dir_to_folder_out,~] = getHALOfileList(site,DATE,'product','cloud');

    % Load, assume only one *_co.nc file per day
    if isempty(halo_files), continue; end
    wstats = load_nc_struct([dir_to_folder_in halo_files{1}]);

    % Create common attribues, fields, and dimensions
    [data,att,dim] = createORcopyCommonAttsDims(wstats,C);
    
    fprintf('\nGenerating the Halo cloud product.\n')

    % Get field names and extract the variables of interest
    fnames = fieldnames(wstats);
    ifs_t = strmatch('time', fnames);
    fnames_time = fnames(ifs_t);

    for ii = 1:length(fnames_time)
        dt_hrs = median(diff((wstats.(fnames_time{ii}))));
        tres = num2str(dt_hrs*60); % time reso in string
        
        fprintf(['HALO cloud product: ' tres ' min resolution - '])
        
        % Weighted data also available?
        weighting = any(ismember(fnames,['radial_velocity_weighted_mean_' tres 'min']));
        
        % Find cloud base height(s)
        cloudbit = get_droplet_bit_matine(wstats.height,wstats.(['beta_mean_' tres 'min']),0);
        A = diff(cloudbit,[],2);
        B = zeros(size(A,1),6); % up to 6 cloud layers
        for k = 1:size(A,1)
            if ~isempty(find(A(k,:)==1))
                lenB = length(find(A(k,:)==1,6));
                B(k,1:lenB) = find(A(k,:)==1,6)+1; % '+1' is for the shift
            end
        end
        
        % Initialize
        cloud_base_velocity = zeros(size(cloudbit,1),6);
        cloud_base_height = zeros(size(cloudbit,1),6);
        if weighting
            cloud_base_weighted_velocity = zeros(size(cloudbit,1),6);
        end
        
        % If cloud free day, leave everything as zero
        if all(cloudbit(:)==0), continue; end
        for i = 1:size(B,1)
            for j = 1:size(B,2)
                if B(i,j) > 0
                    cloud_base_velocity(i,j) = wstats.(['radial_velocity_mean_' tres 'min'])(i,B(i,j));
                    cloud_base_height(i,j) = wstats.height(B(i,j));
                    if weighting
                        cloud_base_weighted_velocity(:,j) = wstats.(['radial_velocity_weighted_mean_' tres 'min'])(i,B(i,j));
                    end
                end
            end
        end
        
        %%-- Add variables --%%
        data.(fnames_time{ii}) = wstats.(fnames_time{ii});
        data.(['cloud_base_velocity_' tres 'min']) = cloud_base_velocity;
        if weighting
            data.(['cloud_base_weighted_velocity_' tres 'min']) = cloud_base_weighted_velocity;
        end
        data.(['cloud_base_height_' tres 'min']) = cloud_base_height;
        data.(['cloud_mask_' tres 'min']) = cloudbit;

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
        % cloud base velocity
        att.(['cloud_base_velocity_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'cloud_layer'},...
            ['Cloud base velocity (' tres ' min)'],...
            {'m s-1','m s<sup>-1</sup>'},...
            C.missing_value,...
            'Discrete time steps',...
            {[-3 3], 'linear'});
        if weighting
            % cloud base weighted velocity
            att.(['cloud_base_weighted_velocity_' tres 'min']) = create_attributes(...
                {['time_' tres 'min'],'cloud_layer'},...
                ['Cloud base weighted velocity (' tres ' min)'],...
                {'m s-1','m s<sup>-1</sup>'},...
                C.missing_value,...
                'Discrete time steps',...
                {[-3 3], 'linear'});
        end
        % cloud base height
        att.(['cloud_base_height_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'cloud_layer'},...
            ['Cloud base height (' tres ' min)'],...
            {'m','m'},...
            C.missing_value,...
            'Discrete time steps',...
            {[0 10000], 'linear'});
        % cloud mask
        att.(['cloud_mask_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'range'},...
            ['Cloud mask (' tres ' min)'],...
            [],...
            C.missing_value,...
            'Simple binary cloud mask, 1: cloud, 0: no cloud.',...
            {[0 1], 'linear'});
        
        % Create dimensions
        dim.(['time_' tres 'min']) = length(data.(['time_' tres 'min']));     
    end
    
    fprintf('\n')
    
    
    % cloud layer
    data.cloud_layer = transpose(1:6);
    att.cloud_layer = create_attributes(...
        {'cloud_layer'},...
        'Cloud layer', ...
        '',...
        [],...
        'Specifies the number of the cloud layer, lowest is the nearest to the instrument. Maximum of six layers can be detected.');
    
    % Add dim
    dim.cloud_layer = size(data.cloud_layer,1);
    
    % Create global attributs
    att.global.Conventions = 'CF-1.0';
    att.global.system = C.system;
    att.global.location = C.location;
    att.global.source = C.source;
    att.global.institution = C.institution;
    att.global.title = C.title;
    att.global.day   = str2num(thedate(7:8));
    att.global.month = str2num(thedate(5:6));
    att.global.year  = str2num(thedate(1:4));
    current_date = datestr(now,'yyyy-mm-dd HH:MM:SS');
    att.global.history = [current_date ' - Created by ' C.user ];
    
    % Order fields
    data = orderfields(data);
    att  = orderfields(att);
    
    % Write into new netcdf
    write_nc_struct(fullfile([dir_to_folder_out '/' thedate '_' site ...
        '_halo-doppler-lidar_cloud.nc']), ...
        dim, data, att)
end



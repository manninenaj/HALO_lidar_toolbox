function calculateHALOwindShearProduct(site,DATES,windproduct,typeof,dt)
%calculateHALOwindShearProduct calculates vector wind shear, and writes 
% the results into a daily *.nc file.
%
% Usage:
% calculateHALOwindShearProduct(site,DATES,windproduct,typeof)
% calculateHALOwindShearProduct(site,DATES,windproduct,typeof,dt)
%
% Inputs:
% -site          String, site name, e.g. site = 'kuopio'
% -DATES         Scalar or vector, numeric, e.g. DATES = 20170401
%                or DATES = [20170401 20170431]
% -windproduct   String, name of the wind product, 'windvad', 'winddbs'
% -elevangle     String, elevation angle 0-90
% -dt            Scalar or vector, numeric, temporal resolution in minutes, 
%                e.g. dt = [1 3 5], or dt = 30, by default dt = [3 5 10 30]
%
% Created 2018-01-18
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

if nargin < 4
  error('''site'', ''DATES'', ''windproduct'', and ''typeof'' are required inputs!')
end
if ~ischar(site)
  error('The 1st input ''site'' must be a string.')
end
if length(DATES)>2
    error('The 2nd input ''DATES'' can have max. length of 2.')
elseif length(DATES)==1
    DATEstart = DATES; DATEend = DATES;
elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
       length(num2str(DATES(2)))~=8)
    error(['The value(s) in the 2nd input ''DATES'' must be' ...
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
    if not(any(strcmp(typeof,{'3','4','5'})))
        error('For dbs winds, the 4th input must be a string and be in the form: ''3'', ''4'', or ''5''.')
    end
end
if nargin < 5
    % Temporal resolutions, min/60 = hrs
    dt = [3 30]; dt = dt./60;
elseif nargin == 5
    if ~isnumeric(dt) | int16(dt)~=dt
        error(['The 5th input must a numerical scalar or vector'...
            ' specifying the temporal resolution in full minutes.'])
    else
        % Temporal resolutions, min/60 = hrs
        dt = dt./60;
    end
else
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
    
    switch windproduct
	case 'windvad'
	    typeof1 = ['ele' typeof];
	case 'winddbs'
            typeof1 = [typeof 'beams'];
    end
    [dir_wind_in, wind_files_tday] = getHALOfileList(site,DATE,'product',windproduct,typeof1);
    if isempty(wind_files_tday)
        switch windproduct
            case 'windvad'
                fprintf('\nNo ''%s-%s'' files found for ''%s'' at ''%s'', skipping...\n',windproduct,typeof1,thedate,site)
            case 'winddbs'
                fprintf('\nNo ''%s-%s'' files found for ''%s'' at ''%s'', skipping...\n',windproduct,typeof1,thedate,site)
        end
        continue
    end
    [~, wind_files_yday] = getHALOfileList(site,DATE_yd,'product',windproduct,typeof1);
    [~, wind_files_tmrw] = getHALOfileList(site,DATE_tw,'product',windproduct,typeof1);

    fprintf('\nGenerating the Halo wind shear product.\n')

    % Get & check output path, is it writable? 
    [dir_wind_shear_out,~] = getHALOfileList(site,DATE,'product','windshear');    
    status = checkHALOpath(site,DATE,'product','windshear');
    if isempty(status)
        fprintf('Cannot write windshear product for the site %s and date %s.',site,num2str(DATE));
        continue;
    end

    % Load, assume only one *.nc file per day, load only the needed fields
    wind_tday = load_nc_struct(fullfile([dir_wind_in '/' wind_files_tday{1}]),{'time','height','u','v','u_error','v_error'});

    % Create common attribues, fields, and dimensions
    [data,att,dim] = createORcopyCommonAttsDims(wind_tday,'product',C);
   
    if ~isempty(wind_files_yday)
        wind_yday = load_nc_struct(fullfile([dir_wind_in '/' wind_files_yday{1}]),{'time','height','u','v','u_error','v_error'});
    else
        % empty otherwise
        wind_yday = [];
    end
    if ~isempty(wind_files_tmrw)
        % load if exists
        wind_tmrw = load_nc_struct(fullfile([dir_wind_in '/' wind_files_tmrw{1}]),{'time','height','u','v','u_error','v_error'});
    else 
        % empty otherwise
        wind_tmrw = [];
    end

    for kk = 1:length(dt)

        fprintf('wind shear: %s min resolution...', num2str(dt(kk).*60))

        % Create time steps and reference
        atime = dt(kk)/2:dt(kk):24 - dt(kk)/2;
        [Xr,Yr] = meshgrid(wind_tday.height, atime);
        
        % Take day before and after to fill some gaps
        if not(isempty(wind_tmrw) || isempty(wind_yday))
            [Xo,Yo] = meshgrid(wind_tday.height,...
                [wind_yday.time(:)-24;...
                wind_tday.time(:);...
                wind_tmrw.time(:)+24]);
            u = interp2(Xo,Yo, ...
                [wind_yday.u; ...
                wind_tday.u; ...
                wind_tmrw.u],Xr,Yr);
            if strcmp(windproduct,'windvad')

            u_error = interp2(Xo,Yo, ...
                [wind_yday.u_error; ...
                wind_tday.u_error; ...
                wind_tmrw.u_error],Xr,Yr);
            end            
            v = interp2(Xo,Yo, ...
                [wind_yday.v; ...
                wind_tday.v; ...
                wind_tmrw.v],Xr,Yr);
            
            if strcmp(windproduct,'windvad')
            v_error = interp2(Xo,Yo, ...
                [wind_yday.v_error; ...
                wind_tday.v_error; ...
                wind_tmrw.v_error],Xr,Yr);
            end
        else
            [Xo,Yo] = meshgrid(wind_tday.height,wind_tday.time);
            u = interp2(Xo,Yo,wind_tday.u,Xr,Yr);
            v = interp2(Xo,Yo,wind_tday.v,Xr,Yr);
            if strcmp(windproduct,'windvad')
            u_error = interp2(Xo,Yo,wind_tday.u_error,Xr,Yr);
            v_error = interp2(Xo,Yo,wind_tday.v_error,Xr,Yr);
            end
        end
        
        % Calculate wind shear
        win_size = 5; % calculate over about 100 meter range
        shear_vec = nan(size(u));
        shear_vec_e = nan(size(u));
        for ir = 1:length(atime)
            for ic = floor(win_size/2)+1:length(wind_tday.height)-floor(win_size/2)
                du = u(ir,ic+floor(win_size/2))-...
                    u(ir,ic-floor(win_size/2));
                dv = v(ir,ic+floor(win_size/2))-...
                    v(ir,ic-floor(win_size/2));

                dz = wind_tday.height(ic + floor(win_size/2))-...
                    wind_tday.height(ic - floor(win_size/2));
                shear_vec(ir,ic) = sqrt((du).^2 + (dv).^2) ./ dz;

                if strcmp(windproduct,'windvad')
                % Error propagation, u_error and v_error are dependent
                shear_vec_e(ir,ic) = sqrt(u_error(ir,ic).^2 +  v_error(ir,ic).^2) ./ dz;
                end
            end
        end
        % time reso in string
        tres = num2str(dt(kk)*60);
        data.(['time_' tres 'min']) = atime(:);
        data.(['vector_wind_shear_' tres 'min']) = shear_vec;
        if strcmp(windproduct,'windvad')
        data.(['vector_wind_shear_error_' tres 'min']) = shear_vec_e;
        end        
        % time
        att.(['time_' tres 'min']) = create_attributes(...
            {['time_' tres 'min']},...
            'Decimal hours UTC', ...
            'Hours UTC',...
            [],...
            ['Discrete time steps, in ' tres ' min temporal resolution.']);
        att.(['time_' tres 'min']).axis = 'T';
        % vector_wind_shear
        att.(['vector_wind_shear_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'range'},...
            'Vector_wind_shear',...
            {'s-1',''},...
            C.missing_value,...
            'Calculated over five range gates which in vertically pointing mode translate to about 100 m.',...
            {[0 0.06], 'linear'});
        if strcmp(windproduct,'windvad')
        % vector_wind_shear error
        att.(['vector_wind_shear_error_' tres 'min']) = create_attributes(...
            {['time_' tres 'min'],'range'},...
            'Vector_wind_shear error',...
            {'s-1',''},...
            C.missing_value,...
            'Propagated from u_error and v_error.',...
            {[0 0.01], 'linear'});
        end
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
    current_date = datestr(now);
    current_date(current_date == '-') = ' ';
    att.global.history = [current_date ' - Created by ' C.user ];
    
    % Order fields
    data = orderfields(data);
    att  = orderfields(att);
    
    % Write into new netcdf
    write_nc_struct(fullfile([dir_wind_shear_out '/' thedate ...
        '_' site '_halo-doppler-lidar_wind-shear.nc']), dim, data, att)
end


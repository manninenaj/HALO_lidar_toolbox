function calculateABLclassificationClimatologyProduct(site,DATES,varargin)
%calculateABLclassificationClimatologyProduct calculates ABL classification 
%climatology product
%
% Usage:
% calculateHALOablClassStats(site,DATES)
%
% Inputs:
% -site          string, site name, e.g. site = 'kuopio'
% -DATES         scalar or vector, numeric, e.g. DATES = 20170401
%                or DATES = [20170401 20170431]
%
% Created 2020-01-14
% Antti Manninen
% Finnish Meteorological Insitute, Finland
% antti.manninen@fmi.fi

% Check inputs
if nargin < 2
    error("'site', and 'DATES' are required inputs!")
end
if ~ischar(site)
    error("The first input 'site' must be a string.")
end
if length(DATES)>2
    error("'DATES' can have max. length of 2.")
elseif length(DATES)==1
    DATEstart = DATES; DATEend = DATES;
elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
        length(num2str(DATES(2)))~=8)
    error(["The value(s) in the second input 'DATES' must be" ...
        " numerical date(s) in YYYYMMDD format."])
else
    DATEstart = DATES(1); DATEend = DATES(2);
end


%%-- Input parameter parsing --%%

% Set defaults
p.height_range_interval = 100; % m agl
p.height_range_max = 4000; % m agl

% Parse inputs
p = parsePropertyValuePairs(p, varargin);

% Validate inputs
if round(p.height_range_max / p.height_range_interval) ~= p.height_range_max / p.height_range_interval
    error("Optional input error: 'height_range_max' has to be divisible with 'height_range_interval'!");
else
    p.num_height_ranges = p.height_range_max / p.height_range_interval;
end

height_bins = 0:p.height_range_interval:p.height_range_max-p.height_range_interval + p.height_range_interval/2;

fprintf('\nGenerating ABL classification climatology product for %s from %i to %i.\n', site, DATEstart, DATEend)
%%-- Load data and calculate --%

% Iterate days
ifirst = true;
fnames_time_out = [];
for daten = datenum(num2str(DATEstart),'yyyymmdd'):...
        datenum(num2str(DATEend),'yyyymmdd')
    
    % Convert date into required formats
    thedate = datestr(daten,'yyyymmdd');
    DATE = str2double(thedate);
    [~,month] = datevec(daten);
    
    % Get default and site/unit/period specific parameters
    C = getconfig(site,DATE);
    
    % Get list of files
    [dir_BLclass,files_BLclass] = getHALOfileList(site,DATE,'product','ABLclassification');
    
    % Check path to write out
    status = checkHALOpath(site,DATE,'level3','ABLclassificationClimatology');
    if isempty(status)
        fprintf('Can''t write %s - %s.',num2str(DATE),site);
        continue;
    end
    [dir_BLc_out,~] = getHALOfileList(site,DATE,'level3','ABLclassificationClimatology');
    
    % Load, assume only one *.nc file per day, and if empty --> skip
    if isempty(files_BLclass), continue; end
    [data, att] = load_nc_struct(fullfile([dir_BLclass, files_BLclass{1}]));
    if isempty(data), continue; end
    
    % Get time resolution(s) from the file
    fnames = fieldnames(data);
    fnames_time = fnames(strmatch('time', fnames));
    
    % ONLY perform at first iteration
    if ifirst
        % extract available time resolution from the the 1st day data,
        % NOTE any other time resolutions found later are collected as well
        fnames_time_out = fnames_time;
        % fetch number of classes from the first time resolution data, same for all..
        tres_1st = fnames_time_out{1}(strfind(fnames_time_out{1},'_')+1:end);
        
        % find number of classes for each product
        str_1 = att.(['bl_classification_' tres_1st]).definition;
        bl_class_types = str2num(transpose(str_1(strfind(str_1, ':')-1)));
        bl_class_types = transpose([transpose(bl_class_types(1)) transpose(bl_class_types(4:end))]);
        bl_class_types(bl_class_types==0) = 2; % shift '0'-no data to '2', '1'-stable/neutral & '2'-unstable omitted
        str_2 = att.(['turbulence_coupling_' tres_1st]).definition;
        turb_cpl_types = str2num(transpose(str_2(strfind(str_2, ':')-1)));
        
        % Initialize output arrays
        for i0 = 1:length(fnames_time_out)
            tres_ini = fnames_time_out{i0}(strfind(fnames_time_out{i0},'_')+1:end);
            % class, month, hour, height
            data_out.(['counts_abl_classification_' tres_ini]) = zeros(numel(bl_class_types), 12, 24, length(height_bins));
            data_out.(['counts_turbulence_coupling_' tres_ini]) = zeros(numel(turb_cpl_types), 12, 24, length(height_bins));
            data_out.(['number_of_elements_' tres_ini]) = zeros(12, 24, length(height_bins));
        end
        
        % set flag for following iterations
        ifirst = false;
        
        % On other iterations check if new time resolutions available
    elseif any(ismember(fnames_time, fnames_time_out))
        fnames_time_new = fnames_time(~ismember(fnames_time, fnames_time_out));
        fnames_time_out = [fnames_time_out; fnames_time_new];
        for ii = 1:length(fnames_time_new)
            tres_new = fnames_time_new{ii}(strfind(fnames_time_out{ii},'_')+1:end);
            % month, hour, range, class
            data_out.(['counts_abl_classification_' tres_new]) = zeros(numel(bl_class_types), 12, 24, length(height_bins));
            data_out.(['counts_turbulence_coupling_' tres_new]) = zeros(numel(turb_cpl_types), 12, 24, length(height_bins));
            data_out.(['number_of_elements_' tres_new]) = zeros(12, 24, length(height_bins));
        end
    end
    
    % Iterate over time resolutions what available for the day in question
    for it = 1:length(fnames_time_out)
        tres = fnames_time_out{it}(strfind(fnames_time_out{it},'_')+1:end);
        
        if ~isfield(data, ['bl_classification_' tres])
            fprintf('\nNo %s resolution data available for %s - %s, skipping.', site, thedate)
            continue
        end
        
        % To select data loop over hrs of day...
        for ihr = 1:24
            % ...and over height ranges
            height_start = 0; height_stop = p.height_range_interval;
            for iheight = 1:length(height_bins)
                
                % Select correct height range and hour
                time_sel = floor(data.(fnames_time{it}))==ihr-1;
                height_sel = data.height_agl>height_start & data.height_agl<=height_stop;
                abl_classification = data.(['bl_classification_' tres])(time_sel, height_sel);
                turbulence_coupling = data.(['turbulence_coupling_' tres])(time_sel, height_sel);
                
                % Remove gaps in data and the 3 lowest most range gates
                %cond = all(abl_classification(:,4:end) == 0, 2);
                %abl_classification(cond,:) = [];
                %abl_classification(:,1:3) = [];
                %turbulence_coupling(cond,:) = [];
                %turbulence_coupling(:,1:3) = [];
                
                % If everything omitted, skip all together.
                if isempty(abl_classification) || isempty(turbulence_coupling)
                    height_start = height_stop;
                    height_stop = height_stop + p.height_range_interval;
                    continue
                end
                
                % shift class '0', i.e. no data, to be '2'
                abl_classification(abl_classification==0) = 2; % 2 = no data from now on
                abl_classification(isnan(abl_classification)) = 2; % 2 = no data from now on
                
                % create edges, abl_classification should be: 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5
                edges_bl_class = [transpose(bl_class_types(:)) bl_class_types(end)+1]-.5;
                
                % calculate counts
                data_out.(['counts_abl_classification_' tres])(:,month,ihr,iheight) = ...
                    squeeze(data_out.(['counts_abl_classification_' tres])(:,month,ihr,iheight)) + ...
                    transpose(histcounts(abl_classification(:), edges_bl_class,'Normalization','count'));
                
                % create edges for turbulence coupling, should be:-0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5
                edges_turb_cpl = [transpose(turb_cpl_types(:)) turb_cpl_types(end)+1]-.5;
                
                % calculate counts
                data_out.(['counts_turbulence_coupling_' tres])(:,month,ihr,iheight) = ...
                    squeeze(data_out.(['counts_turbulence_coupling_' tres])(:,month,ihr,iheight)) + ...
                    transpose(histcounts(turbulence_coupling(:), edges_turb_cpl,'Normalization','count'));
                
                % calculate number of elements
                data_out.(['number_of_elements_' tres])(month,ihr,iheight) = ...
                    squeeze(data_out.(['number_of_elements_' tres])(month,ihr,iheight)) + ...
                    numel(abl_classification(:));
                
                % Update values
                height_start = height_stop;
                height_stop = height_stop + p.height_range_interval;
            end
        end
    end
end


%%-- Add variables --%%
% height
data_out.height_agl = height_bins(:);
att_out.height_agl = create_attributes(...
    {'height_agl'},...
    'Height above ground level', ...
    'm',...
    [],...
    'Center of the height bin.');
att_out.height_agl.standard_name = 'height';

% hours
hr_tmp = .5:1:23.5;
data_out.hour = hr_tmp(:);
att_out.hour = create_attributes(...
    {'hour'},...
    'Hour of the day', ...
    'hr',...
    [],...
    'Center of the bin');
att_out.hour.standard_name = 'hour_of_the_day';

% months
month_tmp = 1:1:12;
data_out.month = month_tmp(:);
att_out.month = create_attributes(...
    {'month'},...
    'Month of the year', ...
    '');
att_out.month.standard_name = 'month_of_the_year';

% bl classification types
data_out.bl_classification_types = bl_class_types(:);
att_out.bl_classification_types = create_attributes(...
    {'bl_classification_types'},...
    'See definition attribute for further details.',...
    '');
att_out.bl_classification_types.definition = ...
    ['2: No signal' 10 ...
    '3: Non-turbulent' 10 ...
    '4: Convective mixing' 10 ...
    '5: Wind shear' 10 ...
    '6: Intermittent' 10 ...
    '7: In cloud' 10 ...
    '8: Cloud-driven' 10 ...
    '9: Precipitation' 10];

% turbulence coupling types
data_out.turbulence_coupling_types = turb_cpl_types(:);
att_out.turbulence_coupling_types = create_attributes(...
    {'turbulence_coupling_types'},...
    'See definition attribute for further details.',...
    '');
att_out.turbulence_coupling_types.definition = ...
    ['0: No signal' 10 ...
    '1: Non-turbulent' 10 ...
    '2: Surface-connected' 10 ...
    '3: Cloud-driven' 10 ...
    '4: In cloud' 10 ...
    '5: Unconnected' 10 ...
    '6: Precipitation' 10];

% Add lat lon
[data_out, att_out] = createLatLonAtts(C, data_out, att_out);

for it2 = 1:length(fnames_time_out)
    tres_out = fnames_time_out{it2}(strfind(fnames_time_out{it2},'_')+1:end);
    %data_out.(['counts_abl_classification_' tres_out]) = histo_data.(['counts_abl_classification_' tres_out]);
    %data_out.(['counts_turbulence_coupling_' tres_out]) = histo_data.(['counts_turbulence_coupling_' tres_out]);
    %data_out.(['number_of_elements_' tres_out]) = histo_data.(['number_of_elements_' tres_out]);
    att_out.(['counts_abl_classification_' tres_out]) = create_attributes(...
        {'bl_classification_types','month','hour','height_agl'},...
        ['Histogram counts calculated from ABL classification ' tres_out ' field per class, month, hour, and every 100 m height range.'],...
        {'',''},...
        C.missing_value);
    att_out.(['counts_turbulence_coupling_' tres_out]) = create_attributes(...
        {'turbulence_coupling_types','month','hour','height_agl'},...
        ['Histogram counts calculated from turbulence coupling ' tres_out ' field per class, month, hour, and every 100 m height range.'],...
        {'',''},...
        C.missing_value);
    att_out.(['number_of_elements_' tres_out]) = create_attributes(...
        {'month','hour','height_agl'},...
        ['Total number of histogram elements in input data (' tres_out ') per month, hour, and every 100 m height range.'],...
        {'',''},...
        C.missing_value);
end

% Create dimensions
dim_out.bl_classification_types = length(data_out.bl_classification_types);
dim_out.turbulence_coupling_types = length(data_out.turbulence_coupling_types);
dim_out.month = length(data_out.month);
dim_out.hour = length(data_out.hour);
dim_out.height_agl = length(data_out.height_agl);

% Create global attributs
att_out.global.Conventions = 'CF-1.0';
att_out.global.system = C.system;
att_out.global.location = C.location;
att_out.global.source = C.source;
att_out.global.institution = C.institution;
att_out.global.title = C.title;
current_date = datestr(now,'yyyy-mm-dd HH:MM:SS');
att_out.global.history = [current_date ' - Created by ' C.user];

% Order fields
data_out = orderfields(data_out);
att_out  = orderfields(att_out);
dim_out  = orderfields(dim_out);

% Write into new netcdf
write_nc_struct(fullfile([dir_BLc_out '/' num2str(DATEstart) '-' num2str(DATEend) '_' site ...
    '_bl-classification_climatology.nc']), dim_out, data_out, att_out)


function [data_out, att_out] = createLatLonAtts(C, data_in, att_in)

att_out = att_in;
data_out = data_in;

% height above...
if isfield(C,'altitude_instrument_level_m_asl') && isfield(C,'altitude_ground_level_m_asl')
    actual_instrument_altitude_asl = C.altitude_instrument_level_m_asl;
    actual_site_altitude_asl = C.altitude_ground_level_m_asl;
    cmnt = '.';
else
    actual_instrument_altitude_asl = C.altitude_in_meters;
    cmnt = ', assumes instrument at ground level becuase no information given.';
    actual_site_altitude_asl = C.altitude_in_meters;
end

% altitude
data_out.altitude_site = actual_site_altitude_asl;
att_out.altitude_site = create_attributes(...
    {},...
    ['Altitude of site above mean sea level' cmnt], ...
    'm');
% altitude of instrument
data_out.altitude_instrument = actual_instrument_altitude_asl;
att_out.altitude_instrument = create_attributes(...
    {},...
    ['Altitude of instrument above mean sea level' cmnt], ...
    'm');

% latitude
if ~isfield(data_in,'latitude')
    if ~isfield(C,'latitude')
        data_out.latitude = 0;
    else
        data_out.latitude = C.latitude;
    end
else
    data_out.latitude = data_in.latitude;
end
att_out.latitude = create_attributes(...
    {},...
    'Latitude of lidar', ...
    'degrees_north');
att_out.latitude.standard_name = 'latitude';

% longitude
if ~isfield(data_in,'longitude')
    if ~isfield(C,'longitude')
        data_out.longitude = 0;
    else
        data_out.longitude = C.longitude;
    end
else
    data_out.longitude = data_in.longitude;
end
att_out.longitude = create_attributes(...
    {},...
    'Longitude of lidar', ...
    'degrees_east');
att_out.longitude.standard_name = 'longitude';

%%


% hf = figure; hf.Units = 'Normalized'; hf.Position = [.2 .1 .4 .65]; hf.Color = 'w';
% % summer
% sp3 = subplot(223);
% bar(N_blclass_summer','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'c)'); title('JJA')
% set(gca,'XTick',2:2:24); xlabel('Time of day UTC'); ylabel('Probability')
% sp3.Position = [.075 .065 .4 .32];
% 
% % autumn
% sp4 = subplot(224);
% bar(N_blclass_autumn','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'd)'); title('SON')
% set(gca,'XTick',2:2:24); xlabel('Time of day UTC');
% sp4.Position = [.55 .065 .4 .32];
% 
% % winter
% sp1 = subplot(221);
% bar(N_blclass_winter','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'a)'); title('DJF')
% set(gca,'XTick',2:2:24); ylabel('Probability')
% sp1.Position = [.075 .48 .4 .32];
% 
% % spring
% sp2 = subplot(222);
% bar(N_blclass_spring','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'b)'); title('MAM')
% set(gca,'XTick',2:2:24);
% sp2.Position = [.55 .48 .4 .32];
% 
% 
% M = {'Missing data','Non-turbulent','Convective mixing','Wind shear',...
%     'Decaying / intermittent','In cloud','Cloud driven'};
% hl = legend(M); hlpos = hl.Position; hlpos(1:2) = [.4 .82];
% hl.Position = hlpos; hl.Box = 'off';
% 
% 
% 
% 
% 
% pause(2)
% export_fig -png -m2 -nocrop blclass_stats_hyytiala_105m_1515m.png
% 
% %%
% 
% hf = figure; hf.Units = 'Normalized'; hf.Position = [.2 .1 .4 .65]; hf.Color = 'w';
% % summer
% sp3 = subplot(223);
% bar(N_epsilon_summer','stacked','EdgeColor','none'); shading flat;
% colormap([epsilon_red(:) epsilon_green(:) epsilon_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'c)'); title('JJA')
% set(gca,'XTick',2:2:24); xlabel('Time of day UTC'); ylabel('Probability')
% sp3.Position = [.075 .065 .4 .32];
% 
% % autumn
% sp4 = subplot(224);
% bar(N_epsilon_autumn','stacked','EdgeColor','none'); shading flat;
% colormap([epsilon_red(:) epsilon_green(:) epsilon_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'd)'); title('SON')
% set(gca,'XTick',2:2:24); xlabel('Time of day UTC');
% sp4.Position = [.55 .065 .4 .32];
% 
% % winter
% sp1 = subplot(221);
% bar(N_epsilon_winter','stacked','EdgeColor','none'); shading flat;
% colormap([epsilon_red(:) epsilon_green(:) epsilon_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'a)'); title('DJF')
% set(gca,'XTick',2:2:24); ylabel('Probability')
% sp1.Position = [.075 .48 .4 .32];
% 
% % spring
% sp2 = subplot(222);
% bar(N_epsilon_spring','stacked','EdgeColor','none'); shading flat;
% colormap([epsilon_red(:) epsilon_green(:) epsilon_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'b)'); title('MAM')
% set(gca,'XTick',2:2:24);
% sp2.Position = [.55 .48 .4 .32];
% 
% 
% M = {'Missing data','Non-turbulent','Connected with surface','Connected with cloud','In cloud','Unconnected'};
% hl = legend(M); hlpos = hl.Position; hlpos(1:2) = [.4 .82];
% hl.Position = hlpos; hl.Box = 'off';
%  
% pause(2)
% export_fig -png -m2 -nocrop epsilon_stats_hyytiala_105m_1515m.png


        % ...hahmotelmia...
        % str_1 = att.(['bl_classification_' dt_i]).definition;
        % bl_class_types = str2num(str_1(strfind(str_1, ':')-1)')
        % bl_class_counts = cell2mat(arrayfun(@(y) sum(ismember(bl_class_tmp,y)')', bl_class_types, 'UniformOutput', 0))
        % str_2 = att.(['turbulence_coupling_' dt_i]).definition;
        % turb_cpl_types = str2num(str_2(strfind(str_2, ':')-1)')
        % turb_cpl_counts = cell2mat(arrayfun(@(y) sum(ismember(turb_cpl_tmp,y)')', turb_cpl_types, 'UniformOutput', 0))
        % counts_reshaped = reshape(counts, p.num_height_ranges, size(counts,2)*p.num_height_ranges)
        % counts_per_hrange = cell2mat(arrayfun(@(x) sum(counts_reshaped(:, x:p.num_height_ranges:end)), [1:p.num_height_ranges], 'UniformOutput', 0)')
        % res = arrayfun(@(x) counts_reshaped(:, x:p.num_height_ranges:end), height_ranges, 'UniformOutput', 0)


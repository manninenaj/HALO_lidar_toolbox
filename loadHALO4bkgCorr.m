function data_4bkgcorr = loadHALO4bkgCorr(site,DATE,measmode)
%loadHALO4bkgCorr loads the required HALO lidar quantities for the 
% background correction and outputs them as a struct variable. Time is 
% converted into decimal hours.
%
% Inputs:
% - site            string, name of the site, e.g. 'kuopio'
% - DATE            scalar, numerical date, e.g. 20171231
% - measmode        string, 'stare', 'ppi', or 'rhi'
% 
% Outputs:
% - data_4bkgcorr   struct
%    .time          vector, time in decimal hours
%    .range         vector, range in original units
%    .snr           array,  uncalibrated signal-to-noise ratio

% version 20171020
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% University of Helsinki, Finland

% Check inputs
if nargin < 3
    error(['''site'', ''DATE'', and ''measmode'' are required inputs!'])
end
if ~ischar(site)
    error('The first input ''site'' must be a string.')
end
if ~isnumeric(DATE) || length(num2str(DATE))~=8 
    error(['The second input ''DATE'' must be a numerical date in' ...
        ' YYYYMMDD format.'])
end
if ~ischar(measmode) || (~strcmp(measmode,'stare') && ...
        ~strcmp(measmode,'ppi') && ~strcmp(measmode,'rhi'))
    error(['The third input ''measmode'' must be a string and can be:'...
         ' ''stare'', ''ppi'', or ''rhi''.'])
end

% Get default and site/unit specific parameters
C = getconfig(site,DATE);

% Convert the date into character array
thedate = num2str(DATE);

% Get file list, always 'uncalibrated' for bkg correction
halo_files = getHALOfileList(site,DATE,'original',measmode);

% Check do data and data fields exist and get field names
fnames_4bgkcorr = checkHALOdata4bkgCorr(site,DATE,measmode);
name_time  = fnames_4bgkcorr{1};
name_range = fnames_4bgkcorr{2};
name_snr   = fnames_4bgkcorr{3};

% initialize
time_tmp = cell(length(halo_files),1);
snr_tmp  = cell(length(halo_files),1);

for i = 1:length(halo_files)
    
    % Load only the time, range, and signal fields
    data = load_nc_struct_silent([C.(['dir_original_' measmode]) ...
        thedate(1:4) '/' halo_files{i}], fnames_4bgkcorr);
    
    % time
    time_tmp{i} = data.(name_time)(:);
    
    % range
    if i==1 % on the first loop
        switch C.range_units
            case 'm'
                range = data.(name_range)(:);
            case 'km'
                range = data.(name_range)(:)/1000;
        end
    end
    
    % snr
    % Check dimensions
    [lrow,lcol] = size(data.(name_snr));
    if lrow == length(data.time)
        snr_tmp{i} = data.(name_snr);
    elseif lcol == length(data.time)
        snr_tmp{i} = data.(name_snr)';
    else
        error(sprintf(['Length of the field ''time''' ...
            ' doesn''t match either dimensions of the' ...
            ' field ''%s'' in\n%s'], name_snr,...
            [directory thedate(1:4) '/' halo_files{i}]))
    end
    
end

% Organize the datasets into one
time_n = cell2mat(time_tmp);
snr = cell2mat(snr_tmp);

% Check format and convert into decimal hours
switch C.time_format
    case 'julian'
        [~,~,~,hh,mm,sss] = jd2date(time_n(:));
        time = hh(:)+mm(:)/60+sss(:)/3600;
    case 'hours'
        time = time_n(:);
    case 'seconds'
        time = time_n(:)/3600;
end

% Outut in a struct
data_4bkgcorr.time = time;
data_4bkgcorr.range = range;
data_4bkgcorr.snr = snr;
end

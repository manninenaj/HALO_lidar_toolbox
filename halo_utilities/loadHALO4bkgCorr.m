function data_4bkgcorr = loadHALO4bkgCorr(site,DATE,measmode,typeof)
%loadHALO4bkgCorr loads the required HALO lidar quantities for the
% background correction and outputs them as a struct variable. Time is
% converted into decimal hours.
%
% Usage:
% data_4bkgcorr = loadHALO4bkgCorr(site,DATE,measmode,typeof)
%
% Inputs:
% - site            string, name of the site, e.g. 'kuopio'
% - DATE            scalar, numerical date, e.g. 20171231
% - measmode        string, 'stare', 'vad', 'dbs', 'rhi', 'custom'
% - typeof          string, 'co','cross','eleXX','aziXXX','type1', '3beams', etc.
%
% Outputs:
% - data_4bkgcorr   struct
%    .time          vector, time in decimal hours
%    .range         vector, range in original units
%    .snr           array,  uncalibrated signal-to-noise ratio
%
% Created 2017-10-20
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% University of Helsinki, Finland

% Check inputs
if nargin < 4
    error(['''site'', ''DATE'', ''measmode'', and ''typeof'' are required inputs!'])
end
if ~ischar(site)
    error('The first input ''site'' must be a string.')
end
if ~isnumeric(DATE) || length(num2str(DATE))~=8
    error(['The second input ''DATE'' must be a numerical date in' ...
        ' YYYYMMDD format.'])
end
if ~ischar(measmode) || ~any(strcmp(measmode,{'stare','vad','dbs','rhi','custom'}))
    error(['The third input ''measmode'' must be a string and can be:'...
        ' ''stare'', ''vad'', ''dbs'', ''rhi'', or ''custom''.'])
end

% Get default and site/unit specific parameters
C = getconfig(site,DATE);

% Get file list, always 'uncalibrated' for bkg correction
[path_to_files, halo_files] = getHALOfileList(site,DATE,'original',measmode,typeof);
if isempty(halo_files)
    data_4bkgcorr = [];
else
    % initialize
    data_4bkgcorr = [];
    time_tmp = cell(length(halo_files),1);
    snr_tmp  = cell(length(halo_files),1);
    
    for i = 1:length(halo_files)
        
        % Load only the time, range, and signal fields
        data = load_nc_struct(fullfile(path_to_files, halo_files{i}), ...
            {C.field_name_original_time,...
            C.field_name_original_range,...
            C.field_name_original_snr});
        if isempty(data)
            error('Can''t find files from %s',...
                [C.(['dir_original_' measmode]) ...
                date_dir])
        end
        
        % Check number of range gates
        if C.num_range_gates ~= length(data.(C.field_name_original_range)(:))
            warning('\n Number of range gates is %s and not %s as specified in the config file --> skipping...',...
                num2str(length(data.(C.field_name_original_range)(:))),num2str(C.num_range_gates))
        else
            % time
            time_tmp{i} = data.(C.field_name_original_time)(:);
            
            % range
            switch C.range_units_original
                case 'm'
                    range = data.(C.field_name_original_range)(:);
                case 'km'
                    range = data.(C.field_name_original_range)(:)/1000;
            end
            
            % snr
            % Check dimensions
            [lrow,lcol] = size(data.(C.field_name_original_snr));
            if lrow == length(data.(C.field_name_original_time))
                snr_tmp{i} = data.(C.field_name_original_snr);
            elseif lcol == length(data.(C.field_name_original_time))
                snr_tmp{i} = transpose(data.(C.field_name_original_snr));
            else
                error(sprintf(['Length of the field ''time''' ...
                    ' doesn''t match either dimensions of the' ...
                    ' field ''%s'' in\n%s'], C.field_name_original_snr,...
                    [directory date_dir halo_files{i}]))
            end
            
            % Organize the datasets into one
            time_n = cell2mat(time_tmp);
            snr = cell2mat(snr_tmp);
            
            % Check format and convert into decimal hours
            switch C.time_format_original
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
    end
end

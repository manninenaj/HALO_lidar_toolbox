function field_names = checkHALOdata4bkgCorr(site,DATE,measmode)
%checkHALOdata4bkgCorr checks that the HALO data file has fields which
% correspond to time, range, signal-to-noise ratio and that their names
% match most commonly used fiels names or a user specified field name.
%
% Inputs:
% - site            string, name of the site, e.g. 'kuopio'
% - DATE            scalar, numerical date, e.g. 20171231
% - measmode        string, 'stare', 'ppi', or 'rhi'
% 
% Outputs:
% - field_names     cell array of strings, field names of time, range, snr
%                   e.g. field_names = {'time'; 'range'; 'snr'}

% version 20171020
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% University of Helsinki, Finland

% Check inputs
if nargin < 3
    error('''site'', ''DATE'', and ''measmode'' are required inputs!')
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

% Get list of files
halo_files = getHALOfileList(site,DATE,'original',measmode);

% Check that all files for the day have the same namaing scheme
for i = 1:length(halo_files)
    
    % Get file info
    fileinfo = nc_info([C.(['dir_original_' measmode]) ...
        thedate(1:4) '/' halo_files{i}]);
    
    % Check if the file was found
    if isempty(fileinfo)
        error('Can''t find files from %s',...
            [C.(['dir_original_' measmode]) thedate(1:4) '/'])
    else
        % Extract field names into a cell array
        [fieldnames{1:length(fileinfo.Dataset),1}] = deal(...
            fileinfo.Dataset.Name);
        % Time
        if strcmp(C.fieldname_time,'-') % if not specified for the site
            time_fields = {'time'};
        else 
            time_fields = C.fieldname_time;
        end
        if ~any(ismember(fieldnames,time_fields))
            error(sprintf(['No field which would correpond to ''time'''...
                ' in \n%s'],fileinfo.Filename))
        else
            % identify the time field
            itimefn = ismember(fieldnames,time_fields);
            % take field name from the first file, if more than one
            if i == 1, time_fieldname = fieldnames{itimefn}; end
        end
        % Range
        if strcmp(C.fieldname_range,'-') % if not specified for the site
            range_fields = {'range'};
        else 
            range_fields = C.fieldname_range;
        end
        if ~any(ismember(fieldnames,range_fields))
            error(sprintf(['No field which would correspond to' ...
                ' ''range from the instruments'' and whose name would' ...
                ' match with one of these:\n' repmat('''%s''\n',1,...
                length(range_fields)) 'in %s'], range_fields{:}, ...
                fileinfo.Filename))
        else
            % identify the range field
            irangefn = ismember(fieldnames,range_fields);
            % take field name from the first file, if more than one
            if i == 1, range_fieldname = fieldnames{irangefn}; end
        end
        % Signal-to-noise ratio
        if strcmp(C.fieldname_SNR,'-') % if not specified for the site
            snr_fields = {'snr','SNR','signal','intensity'};
        else
            snr_fields = C.fieldname_SNR;
        end
        if ~any(ismember(fieldnames,snr_fields))
            error(sprintf(['Can''t find a field which would' ...
                ' correspond to ''signal-to-noise ratio'' and whose'...
                ' name would match with one of these:\n' repmat(...
                '''%s''\n',1,length(snr_fields)) 'in %s'],...
                snr_fields{:},fileinfo.Filename))
        elseif sum(ismember(fieldnames,snr_fields)) > 1
            % If field name for signal-to-noise ratio is not given in the
            % config file, more than one of the 'guessed' names can occur.
            % In this case ask for the exact field name from the user.
            error(sprintf(['There is more than one data field whose name'...
                ' can be associated with ''signal-to-noise ratio''.\n'...
                'Give the exact field name in the halo_config.txt'...
                ' file as a parameter ''fieldname_SNR''.\n'...
                'By default these are looked for:\n' repmat(...
                '''%s''\n',1,length(snr_fields))],snr_fields{:}))
        else
            % identify the signal-to-noise ratio field
            isnrfn = find(ismember(fieldnames,snr_fields),1,'last');
            % take field name from the first file, if more than one
            if i == 1, snr_fieldname = fieldnames{isnrfn}; end
        end
    end
    % Put field names into a cell array
    field_names = {time_fieldname;range_fieldname;snr_fieldname};
end

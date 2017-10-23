function file_list = getHALOfileList(site,DATE,processlev,measmode)
%getHALOfileList generates a list of HALO file names measured with a given
% measurement mode and with a given processing level. 
%
% Inputs:
% - site            string, name of the site, e.g. 'kuopio'
% - DATE            scalar, numerical date, e.g. 20171231
% - measmode        string, 'stare', 'ppi', 'rhi', 'vertical', or 'wind'
% - processlev      string, 'uncalibrated', or 'calibrated'
% 
% Outputs:
% - file_list       cell array of string,  file names for the site, 
%                   measurement mode and processing level

% version 20171020
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% University of Helsinki, Finland

% Check inputs
if nargin < 4
    error(['''site'', ''DATE'', ''processlev'', and ''measmode'' are' ...
        ' required inputs!'])
end
if ~ischar(site)
    error('The first input ''site'' must be a string.')
end
if ~isnumeric(DATE) || length(num2str(DATE))~=8
    error(['The second input ''DATE'' must be a numerical date in' ...
        ' YYYYMMDD format.'])
end
if ~ischar(processlev) || (~strcmp(processlev,'original') && ...
        ~strcmp(processlev,'corrected') &&  ...
        ~strcmp(processlev,'calibrated'))
    error(['The third input ''processlev'' must be a string and can be:'...
        ' ''original'', ''corrected'', or ''calibrated''.'])
end
if ~ischar(measmode) || (~strcmp(measmode,'stare') && ...
        ~strcmp(measmode,'ppi') && ~strcmp(measmode,'rhi') && ...
        ~strcmp(measmode,'vertical') && ~strcmp(measmode,'wind'))
    error(['The third input ''measmode'' must be a string and can be:'...
         ' ''stare'', ''ppi'', ''rhi'' ,''vertical'', or ''wind''.'])
end

% Get default and site/unit specific parameters
C = getconfig(site,DATE);

% Initialize
file_list = {};

% Convert the date into character array
thedate = num2str(DATE);

% check if path for given combination of 'processlev' and 'measmode' exist
cpmfield = ['dir_' processlev '_' measmode]; % C.processlev_measmode field
if ~isfield(C,cpmfield)
    error(['Can''t find parameter ''%s'' for the site ''%s'' \nand'...
        ' which would be valid for the date ''%s'' from halo_config.txt'...
        ' file.'], cpmfield,site,num2str(DATE))
end

% Generate path
switch processlev
    case {'original','corrected'}
        % correction keeps the file format and naming the same
        fileformatfield = C.(['file_format_original_' measmode]);
        filenamingfield = C.(['file_naming_original_' measmode]);
    otherwise
        % Assume calibrated files follow Cloundet naming scheme
        % and have format *.nc
        fileformatfield = '.nc';
end
        
dir_to_folder = [C.(cpmfield), thedate(1:4) '/'];
file_names_2look4 = ['*' thedate '*' fileformatfield];
direc = dir([dir_to_folder,file_names_2look4]);
if isempty(direc)
    error(sprintf(['No *%s files for the date %s in\n%s\nCheck paths' ...
        ' in halo_config.txt file.'],fileformatfield,...
        thedate,dir_to_folder));
else
    % Get list of files
    [file_list{1:length(direc),1}] = deal(direc.name);
    switch processlev
        case {'original','corrected'}
            % Look for specific files based on the naming 
            i_file = ~cellfun('isempty',strfind(file_list,...
                filenamingfield));
            file_list = file_list(i_file);
            file_list = sort(file_list);
        otherwise
            % If calibrated, assume Cloudnet naming schmeme and that there 
            % are no other measurement in the same folder
            file_list = sort(file_list);
    end
end
end

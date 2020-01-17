function [dir_to_folder, file_list] = getHALOfileList(site,DATE,processing_level,observation_type,sub_type)
%getHALOfileList generates a list of HALO file names measured with a given
% measurement mode and with a given processing level. Also outputs full path to files.
%
% Usage:
% [~, file_list] = getHALOfilelist(site,DATE,processing_level,observation_type)
% [dir_to_folder,~] = getHALOfileList(site,DATE,processing_level,observation_type)
% [dir_to_folder, file_list] = getHALOfileList(site,DATE,processing_level,observation_type)
% [dir_to_folder, file_list] = getHALOfileList(site,DATE,processing_level,observation_type,sub_type)
%
% Inputs:
% - site              string, name of the site, e.g. 'kuopio'
% - DATE              scalar, numerical date, e.g. 20171231
% - processing_level  string, e.g. 'calibrated'
% - observation_type  string, e.g. 'stare', 'epsilon'
% - sub_type          string, e.g. 'co', 'cross', '3beams', 'ele15', 'co12', 'azi270')      
%
% Outputs:
% - dir_to_folder   full path to the folder 
% - fille_list      cell array of string,  file names for the site,
%                   measurement mode and processing level
%
% Created 2017-10-20, updated 2020-01-17
% Antti Manninen
% antti.manninen@hfmi.fi
% Finnish Meteorological Institute


% Check inputs
list_of_processing_levels = {'original','corrected','calibrated','background','product','level3'};
list_of_observation_types = {'stare','vad','dbs','rhi','custom','co','txt','nc'}; % TODO: make txt and nc optional inputs
list_of_products = {'windvad','winddbs','epsilon','wstats','wstats4precipfilter','sigma2vad','windshear',...
    'LLJ','ABLclassification','cloud','betavelocovariance','ABLclassificationClimatology'};
if nargin < 4
    error("At least inputs 'site', 'DATE', 'processing_level', and 'observation_type'")
end
if (nargin == 4 || nargin == 5) && (strcmp(processing_level,'product') && any(strcmp(observation_type,list_of_products)) || ...
       strcmp(processing_level,'background'))
    if ~ischar(site)
        error("The 1st input (site name) must be a string.")
    end
    if ~isnumeric(DATE) || length(num2str(DATE))~=8
        error("The 2nd input (date) must be a numeric value in YYYYMMDD format.")
    end
    if ~ischar(processing_level) || ~any(strcmp(processing_level, list_of_processing_levels))
        error("The 3rd input (processing level) must be a string and one of these:\n%s", ...
           sprintf('%s,', list_of_processing_levels{:}))
    end
    if ~ischar(observation_type) || ~any(strcmp(observation_type,[list_of_observation_types, list_of_products]))
        error("The 4th input (observation type) must be a string and one of these:\n%s", ...
            sprintf("'%s','%s',", list_of_observation_types{:}, list_of_products{:}))
    end
end
if nargin < 5 && (~strcmp(processing_level,'product') && ~any(strcmp(observation_type, list_of_products)) && ...
       ~strcmp(processing_level,'background'))
    error("Input 'sub_type' is required with %s'", observation_type)
end

% Get default and site/unit specific parameters
C = getconfig(site,DATE);

% Initialize
file_list = {};

% Convert the date into character array
thedate = num2str(DATE);

% check if path for given combination of 'processing_level' and 'observation_type' exist
switch nargin
  case 5
    cpmt = ['dir_' processing_level '_' observation_type '_' sub_type];
  case 4
    cpmt = ['dir_' processing_level '_' observation_type];
end
if ~isfield(C,cpmt)
    error(["Can't find parameter '%s' for the site '%s'\n"...
        "that would be valid on '%s', as specified in halo_config.txt."], ...
        cpmt,site,num2str(DATE))
end

% Get directory to the folder
dir_to_folder = C.(cpmt);

% Replace/expand dates subdirectories as required
% find year, month, day wildcards
dir_to_folder = strrep(dir_to_folder,'+YYYY+',thedate(1:4));
dir_to_folder = strrep(dir_to_folder,'+MM+',thedate(5:6));
if not(isempty(strfind(dir_to_folder,'+DOY+')))
    dir_to_folder = strrep(dir_to_folder,'+DOY+', my_doy(thedate));
else
    dir_to_folder = strrep(dir_to_folder,'+DD+',thedate(7:8));
end

% Generate path
switch processing_level
    case 'original'
     file_naming = C.(['file_naming_original_' observation_type '_' sub_type]);
     file_format = C.file_format_after_hpl2netcdf;
     %   switch ~isempty(findstr(site,'arm-'))
     %    case 1 % ARM site
     %     file_names_2look4 = ['*' file_naming '*' thedate '*' file_format];
     %    case 0 % non-ARM site
     file_names_2look4 = ['*' thedate '*' file_naming  '*' file_format];
     % Check if empty , try other pattern     
     if isempty(dir([dir_to_folder,file_names_2look4]))
         file_names_2look4 = ['*' file_naming '*' thedate '*' file_format];
     end
     %    end
  case 'background'
    switch observation_type
      case 'txt'
        file_format = '.txt';
        file_names_2look4 = ['*' thedate([7:8 5:6 3:4]) '*' ...
                            file_format];
      case 'nc'
        file_format = '.nc';
        file_names_2look4 = ['*' thedate '*' file_format];
    end
 otherwise
  % blindly assume that there are no other type of files
  % with the same naming in the same directory
  file_format = '.nc';
  file_names_2look4 = ['*' thedate '*' file_format];
end

direc = dir([dir_to_folder,file_names_2look4]);
if isempty(direc)
  file_list = [];
  if nargin < 5 && nargout == 1
      warning(["Can't find " processing_level " " observation_type " files for site " site " and date " num2str(DATE) "!"])
  elseif nargout == 1
      warning(["Can't find " processing_level " " observation_type " " sub_type " files for site " site " and date " num2str(DATE) "!"])
  end      
else
    % Get list of files
    [file_list{1:length(direc),1}] = deal(direc.name);
    switch processing_level
        case 'original'
            % More complex for AMR naming scheme..
            if length(site)>3 && strcmp(site(1:4),'arm-')
                % assume ARM file naming scheme, find data level indicator from the name
                b = cellfun(@(x) strfind(x,'.'),file_list,'UniformOutput',false);
                b1sts = cellfun(@(x) x(1),b)+2; % location of 1st dot plus 2 is data level
                fnids = nan(length(b1sts),1); % initialize
                for i = 1:length(b1sts),fnids(i) = str2double(file_list{i}(b1sts(i))); end
                % Select the files with lowest data level, i.e. original
                theuniqs = unique(fnids);
                if numel(theuniqs)>1
                    [~,imin] = min(theuniqs);
                    file_list = file_list(fnids == theuniqs(imin));
                end
                i_file = ~cellfun('isempty',strfind(file_list,...
                    C.(['file_naming_original_' observation_type '_' sub_type])));
                file_list = file_list(i_file);
                file_list = sort(file_list);
            else
                %Look for specific files based on the naming
                i_file = ~cellfun('isempty',strfind(file_list,...
                    file_naming));
                file_list = file_list(i_file);
                file_list = sort(file_list);
            end
        case 'corrected'
            % More complex for AMR naming scheme..
            if length(site)>3 && strcmp(site(1:4),'arm-')
                % assume ARM file naming scheme, find data level indicator from the name
                b = cellfun(@(x) strfind(x,'.'),file_list,'UniformOutput',false);
                b1sts = cellfun(@(x) x(1),b)+2; % location of 1st dot plus 2 is data level
                fnids = nan(length(b1sts),1); % initialize
                for i = 1:length(b1sts), fnids(i) = str2double(file_list{i}(b1sts(i))); end
                % Select the files with highest data level, i.e. latest
                theuniqs = unique(fnids);
                if numel(theuniqs)>1
                    [~,imax] = max(theuniqs);
                    file_list = file_list(fnids == theuniqs(imax));
                end
                i_file = ~cellfun('isempty',strfind(file_list,...
                    C.(['file_naming_original_' observation_type '_' sub_type])));
                file_list = file_list(i_file);
                file_list = sort(file_list);
            else
                %Look for specific files based on the naming
                i_file = ~cellfun('isempty',strfind(file_list,...
                    file_naming));
                file_list = file_list(i_file);
                file_list = sort(file_list);
            end
        otherwise
            % blindly assume that there are no other type of files
            % with the same naming in the same directory
            file_list = sort(file_list);
    end
end
end

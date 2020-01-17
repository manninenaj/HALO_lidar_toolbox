function status = checkHALOpath(site,DATE,processing_level,observation_type,sub_type)
%checkHALOpath checks do the date subdirectories exist, if not create
%if they are located at the end of the stem path.
%
% Usage:
% status = checkHALOpath(site,DATE,processlev,measmode)
%
% Inputs:
% - site              string, name of the site, e.g. 'kuopio'
% - DATE              scalar, numerical date, e.g. 20171231
% - processing_level  string, e.g. 'calibrated'
% - observation_type  string, e.g. 'stare', 'epsilon'
% - sub_type          string, e.g. 'co', 'cross', '3beams', 'ele15', 'co12', 'azi270')
%
% Outputs:
% - status      1 if path exits or was created succesfully, empty otherwise
%
% Created 2017-11-23, updated 2020-01-17
% Antti Manninen
% antti.manninen@fmi.fi
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


% Get default and site specific parameters
C = getconfig(site,DATE);

% Convert to string
thedate = num2str(DATE);

% check if path for given combination of 'processing_level' and 'observation_type' exist
if nargin == 4
    cpmt = ['dir_' processing_level '_' observation_type]; % -C-.-p-rocesslev_-m-easmode
elseif nargin == 5
    cpmt = ['dir_' processing_level '_' observation_type '_' sub_type]; % -C-.-p-rocesslev_-m-easmode -t-ypeof
end
if ~isfield(C,cpmt)
    error(['Can''t find parameter ''%s'' for the site ''%s'' \nand'...
        ' which would be valid for the date ''%s'' from halo_config.txt'...
        ' file.'], cpmt,site,num2str(DATE))
end

% Replace/expand dates subdirectories as required
% find year, month, day wildcards
output_path = C.(cpmt);
dir_to_folder = output_path;
dir_to_folder = strrep(dir_to_folder,'+YYYY+',thedate(1:4));
dir_to_folder = strrep(dir_to_folder,'+MM+',thedate(5:6));
dir_to_folder = strrep(dir_to_folder,'+DD+',thedate(7:8));

if exist(dir_to_folder,'dir')
  status = 1; 
  return;
end

% Strip date subdirectories
pathstem = output_path;
pathstem = strrep(pathstem,'+YYYY+','');
pathstem = strrep(pathstem,'+MM+','');
pathstem = strrep(pathstem,'+DD+','');
pathstem = fullfile(pathstem);

% Simple case (all date subdirectories at end of path)
if strcmp(pathstem,output_path(1:length(pathstem)))
  [status_logic, ~, ~] = mkdir(dir_to_folder);
  if status_logic
    status = 1;
  else
    status = 0;
  end
else
  status = 1;
  fprintf('No directories created.\n');
end

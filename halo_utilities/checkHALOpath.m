function status = checkHALOpath(site,DATE,processlev,measmode,typeof)
%checkHALOpath checks do the date subdirectories exist, if not create
%if they are located at the end of the stem path.
%
% Usage:
% status = checkHALOpath(site,DATE,processlev,measmode)
%
% Inputs:
% - site        string, name of the site, e.g. 'kuopio'
% - DATE        scalar, numerical date, e.g. 20171231
% - processlev  string, 'corrected','calibrated','product'
% - measmode    string, 'stare','ppi','rhi','co','windvad','winddbs',
%               'velostats','epsilon','sigma2vad','BLclassification','cloudmask'
%
% Outputs:
% - status      1 if path exits or was created succesfully, empty otherwise
%
% Created 2017-11-23
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% University of Helsinki, Finland

% Check inputs
if nargin < 5 && (~strcmp(processlev,'product') && ~any(strcmp(measmode,{'epsilon',...
        'wstats','wstats4precipfilter','sigma2vad','windshear','LLJ','ABLclassification','cloud','betavelocovariance'})))
    error(sprintf(['Inputs ''site'', ''DATE'', ''processlev'', ''measmode'', and ''typeof'''...
        ' are required for any other products than: \n''epsilon'', ''wstats'', ''wstats4precipfilter''', ...
        ' ''sigma2vad'',''windshear'', ''LLJ'', ''ABLclassification'', ''cloud'',''betavelocovariance''']))
    if ~ischar(site)
        error('The 1st input ''site'' must be a string.')
    end
    if ~isnumeric(DATE) || length(num2str(DATE))~=8
        error(['The 2nd input ''DATE'' must be a numerical date in' ...
            ' YYYYMMDD format.'])
    end
    if ~ischar(processlev) || ~any(strcmp(processlev,{'original','corrected','calibrated','background','product'}))
        error(['The 3rd input ''processlev'' must be a string and can be:'...
            ' ''original'', ''corrected'', ''calibrated'', ''background'', or ''product''.'])
    end
    if ~ischar(measmode) || ~any(strcmp(measmode,{'stare','vad','dbs','rhi','co','windvad','winddbs',...
            'txt','wstats','wstats4precipfilter','epsilon','sigma2vad','windshear','LLJ','ABLclassification','cloud','betavelocovariance'}))
        error(sprintf(['The 4th input ''measmode'' must be a string and can be:\n'...
            '''stare'',''vad'',''rhi'',''dbs'',''co'',''windvad'',''winddbs'',''txt'',''wstats''\n'...
            '''wstats4precipfilter'', ''epsilon'',''sigma2vad'',''windshear'',''LLJ'',''ABLclassification'',''cloud'',''betavelocovariance''.']))
    end
else
    if nargin < 4
        error(sprintf(['Inputs ''site'', ''DATE'', ''processlev'', ''measmode'''...
            ' are required for the products: \n''epsilon'', ''wstats'', ''wstats4precipfilter'', ''sigma2vad''',...
            '''windshear'', ''LLJ'', ''ABLclassification'', ''cloud'',''betavelocovariance''']))
    end
    if ~ischar(site)
        error('The 1st input ''site'' must be a string.')
    end
    if ~isnumeric(DATE) || length(num2str(DATE))~=8
        error(['The 2nd input ''DATE'' must be a numerical date in' ...
            ' YYYYMMDD format.'])
    end
    if ~ischar(processlev) || ~any(strcmp(processlev,{'original','corrected','calibrated','background','product'}))
        error(['The 3rd input ''processlev'' must be a string and can be:'...
            ' ''original'', ''corrected'', ''calibrated'', ''background'', or ''product''.'])
    end
    if ~ischar(measmode) || ~any(strcmp(measmode,{'stare','vad','dbs','rhi','co','custom','windvad','winddbs',...
            'txt','wstats','wstats4precipfilter','epsilon','sigma2vad','windshear','LLJ','ABLclassification','cloud','betavelocovariance'}))
        error(sprintf(['The 4th input ''measmode'' must be a string and can be:\n'...
            '''stare'',''vad'',''dbs'',''rhi'',''co'',''custom'',''windvad'',''winddbs'',''txt'',''wstats''\n'...
            '''wstats4precipfilter'',''epsilon'',''sigma2vad'',''windshear'',''LLJ'',''ABLclassification'',''cloud'',''betavelocovariance''.']))
    end
end

% Get default and site specific parameters
C = getconfig(site,DATE);

% Convert to string
thedate = num2str(DATE);

% check if path for given combination of 'processlev' and 'measmode' exist
if nargin == 4
    cpmt = ['dir_' processlev '_' measmode]; % -C-.-p-rocesslev_-m-easmode
elseif nargin == 5
    cpmt = ['dir_' processlev '_' measmode '_' typeof]; % -C-.-p-rocesslev_-m-easmode -t-ypeof
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

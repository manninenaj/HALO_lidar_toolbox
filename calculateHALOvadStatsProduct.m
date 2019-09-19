function calculateHALOvadStatsProduct(site,DATES,varargin)
%calculateHALOvadStatsProduct
%
% Usage:
% calculateHALOvadStatsProduct(site,DATES)
% calculateHALOvadStatsProduct(site,DATES,'dt',60)
% calculateHALOvadStatsProduct(site,DATES,'dt',[60 120 240])
%
% Created 2019-09-10
% Antti Manninen
% Finnish Meteorological Institute
% antti.manninen@fmi.fi

if nargin < 3
  error('''site'', ''DATES'', ''elevation'' are required inputs!')
end
  if ~ischar(site)
  error('The first input ''site'' must be a string.')
end
  if length(DATES)>2
  error('''DATES'' can have max. length of 2.')
  elseif length(DATES)==1
  DATEstart = DATES; DATEend = DATES;
elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && length(num2str(DATES(2)))~=8)
  error(['The value(s) in the second input ''DATES'' must be numerical date(s) in YYYYMMDD format.'])
else
  DATEstart = DATES(1); DATEend = DATES(2);
end
if (~ischar(elevangle) || length(elevangle) ~= 2 || (~isempty(str2num(elevangle)) && str2num(elevangle)<0 || str2num(elevangle)>90)) & not(strcmp(elevangle,'0'))
  error('The 3rd input must be a string and no longer than 2 characters specifying the elevation angle 0-90 degrees.')
end

for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):datenum(num2str(DATEend),'yyyymmdd')

  % Convert date into required formats
  thedate = datestr(DATEi,'yyyymmdd');
  DATE = str2double(thedate);

  % Get default and site/unit/period specific parameters
  C = getconfig(site,DATE);

  % Get list of files
  elevangle1 = ['ele' elevangle];
  abc = ['vad_' elevangle1];
  [dir_vad,files_vad] = getHALOfileList(site,DATE,'calibrated','vad',elevangle1);
  [dir_out,~] = getHALOfileList(site,DATE,'product','vadstats',elevangle1);
  [dir_wind,files_wind] = getHALOfileList(site,DATE,'product','windvad',elevangle1);
 
  if isempty(files_vad)
     continue;
  end

  % Check path to write out
  status = checkHALOpath(site,DATE,'product','windvad',elevangle1);
  if isempty(status)
    fprintf('Can''t write %s - %s.',num2str(DATE),site);
    continue;
  end

  % Read winds, one file per day, one netcdf file
  [datawind,~,~] = load_nc_struct(fullfile([dir_wind '/' files_wind{i}]));

  % Initialise
  signal = cell(length(files_vad));

  % One file per scan, several scans per day, many netcdf files
  for i = 1:length(files_vad)
    % Load
    [datavad,~,~] = load_nc_struct(fullfile([dir_vad '/' files_vad{i}]));

    % Quick & dirty clutter-noise map, important especially in urban environments
    % Better, average several scans, noise should average to "zero", 
    % strong returns from hard targets should become obvious and follow gaussian,
    % find hard targer location via convolution and gaussian filter   
    myfilter = datavad.signal > 1.2 | datavad.signal < 1.01;  
    datavad.signal(myfilter) = nan;
    datavad.beta_raw(myfilter) = nan;
    datavad.v_raw(myfilter) = nan;

    % Normalize radial velocity with mean wind
    v_norm = datavad.v_raw ./ repmat(datawind.wind_speed(i,:),length(datavad.azimuth),1);
    
    % Collect to cells  
    signal_cell{i} = permute(datavad.signal, [3 1 2]);
    beta_cell{i} = permute(datavad.beta_raw, [3 1 2]);
    v_raw_cell{i} = permute(datavad.v_raw, [3 1 2]);
    v_error_cell{i} = permute(datavad.v_error, [3 1 2]);)
    v_norm_cell{i} = permute(v_norm, [3 1 2]);
    time_cell{i} = nanmedian(decimal2daten(datavad.time,DATEi));
    azim_cell{i} = datavad.azimuth; 
  end
    
  % dims: scantime x azimuth x range
  signal = cell2mat(signal_cell(:));
  beta = cell2mat(beta_cell(:)); 
  velo = cell2mat(v_raw_cell(:));
  velo_error = cell2mat(v_error_cell(:));
  velo_norm = cell2mat(v_norm_cell(:));
  time = cell2mat(time_cell);
  azim = nanmedian(cell2mat(azim_cell));
  
  % average and with mean wind
  for idt = 1:length(dt)
    starttime = DATEi;
    endtime = DATEi + datenum(0,0,0,0,dt(idt),0);
    inds = find(time>=starttime && time<=endtime);
    
    velo_mean = nanmean(squeeze(velo(inds,:,:)0);
    velo_var = nanvar(squeeze(velo(inds,:,:)));
    velo_var_unbiased = nanvar(squeeze(velo(inds,:,:))) - nanvar(squeeze(velo_error(inds,:,:)));










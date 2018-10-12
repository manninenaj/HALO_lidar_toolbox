function extractHALOwstatsNormalization(site,DATES)
%
%
%
% Usage:
% - one day only
% extractHALOwstatsNormalization('kuopio',20160827)
% - from to
% extractHALOwstatsNormalization('kuopio',[20160827 20160830])
% - selected dates
% extractHALOwstatsNormalization('kuopio',[20160827 20160828 20160830 20161011])
%
%
% Check inputs
DATEstart = []; DATEend = []; DATES_cell = []; 
if nargin < 2
    error('''site'' and ''DATES'' are required inputs!')
end
if ~ischar(site) && ~isempty(site)
    error('The first input ''site'' must be a string, and can''t be empty.')
end
if ~isempty(DATES)
    if length(DATES)==1
        iter_flag = 1;
        if length(num2str(DATES))~=8
            error(['The value in the second input ''DATES'' must be' ...
                ' numerical date in YYYYMMDD format.'])
        end
    elseif length(DATES)==2
        iter_flag = 2;
        DATEstart = DATES(1); DATEend = DATES(2);
        if length(num2str(DATEstart))~=8 || length(num2str(DATEend))~=8
            error(['The second input ''DATES'' must be' ...
                ' numerical date(s) in YYYYMMDD format, e.g. [20180627 20180630]'])
        end
    else
        iter_flag = 3;
        DATES_cell = cellstr(num2str(DATES(:)));
        if ~all(cellfun(@(x) length(x)==8, DATES_cell))
            error(['The second input ''DATES'' must be' ...
                ' numerical date(s) in YYYYMMDD format.'])
        end
    end
else
    error('The second input ''DATES'' can''t ne empty.')
end

switch iter_flag
    case 1
        DATES_iter = datenum(num2str(DATES),'yyyymmdd');
    case 2
        DATES_iter = datenum(num2str(DATEstart),'yyyymmdd'):...
            datenum(num2str(DATEend),'yyyymmdd');
    case 3
        DATES_iter = datenum(DATES_cell,'yyyymmdd');
end


bkg_var_cell = cell(length(DATES_iter),1);
for i = 1:length(DATES_iter)
    % Check input and ouput files
    thedate = datestr(DATES_iter(i),'yyyymmdd');
    DATE = str2double(thedate);
        
    % Get list of wstats files
    [dir_in, files] = getHALOfileList(site,DATE,'product','wstats4precipfilter');
    % If no files for today, skip the day
    if isempty(files)
        fprintf('\nNo ''wstats4precipfilter'' files found for ''%s'' at ''%s'', skipping...\n',thedate,site)
        continue
    end
    
    % Get & check output path, can the file be written
    [dir_out,~] = getHALOfileList(site,DATE,'product','wstats4precipfilter');
    status = checkHALOpath(site,DATE,'product','ABLclassification');
    if isempty(status)
        fprintf('Cannot write wstats4precipfilter normalization file.');
        continue;
    end

    datab = load_nc_struct(fullfile([dir_in files{1}]),...
        {'beta_variance_3min','beta_wvariance_3min','height'});
    
    % Collect 
    bkg_var_cell{i} = datab.beta_variance_3min;

end

% Convert from cell to mat
bkg_var = cell2mat(bkg_var_cell);

% Calculate mean bkg variance, assume it follows normal distribution
bkg_var_mean = nanmedian(bkg_var);

% Initialize
Output = cell(1,2);

% Construct a cell array which will be written to output file
for iCellRow = 1:length(datab.height)
    Output{iCellRow,1} = datab.height(iCellRow);
    Output{iCellRow,2} = bkg_var_mean(iCellRow) / 1e-13;
    Output{iCellRow,3} = '1e-13';
end

% Open file for writing
if length(DATES_iter)>1
    fileID = fopen(sprintf('%s%s_%s_background_mean_variance_profile.dat',...
        dir_out, datestr(DATES_iter(1),'yyyymmdd'), datestr(DATES_iter(2),'yyyymmdd')), 'w');
else
    fileID = fopen(sprintf('%s%s_%s_background_mean_variance_profile.dat',...
        dir_out, datestr(DATES_iter(1),'yyyymmdd'), datestr(DATES_iter(2),'yyyymmdd')), 'w');
end

% Specify format
formatSpec = '%.13g %.13g %s\n';

% Write to file row by row
for iRow = 1:length(datab.height), fprintf(fileID, formatSpec, Output{iRow,:}); end

% Close file
fclose(fileID);

end



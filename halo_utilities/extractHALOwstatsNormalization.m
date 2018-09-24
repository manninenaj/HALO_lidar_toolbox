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
DATEstart = []; DATEend = []; DATES_cell = []; iter_flag = [];
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

for DATEi = DATES_iter
    % Check input and ouput files
    thedate = datestr(DATEi,'yyyymmdd');
    DATE = str2double(thedate);
    
    % Get default and site/unit/period specific parameters
    C = getconfig(site,DATE);
    
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

    % load signal and beta
    data = load_nc_struct(fullfile([dir_in '/' files{1}]));

    disp('Hello!')
end



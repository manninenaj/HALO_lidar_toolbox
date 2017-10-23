function runHALObackgroundCorrection(site,DATES)
%runHALObackgroundCorrection loads vertically pointing (stare) and plan 
% position indicator (PPI) HALO data, corrects the background signal, and 
% writes the corrected data into new separate stare and PPI files while 
% keeping the original number of files per day and the original data field 
% structure of the files.
%
% Inputs:
% - site         string, name of the site, e.g. site = 'sodankyla'
% - DATES        scalar or vector, numerical date in decimal hours, 
%                e.g. DATES = 20170101 or DATES = [20170101,20170131]

% version 20171019
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% Department of Physics
% University of Helsinki, Finland

% Check inputs
if nargin < 2
    error('''site'' and ''DATES'' are required inputs!')
end
if ~ischar(site)
    error('The first input ''site'' must be a string.')
end
if length(DATES)>2
    error('''DATES'' can have max. length of 2.')
elseif length(DATES)==1
    if length(num2str(DATES))~=8
    error(['The value in the second input ''DATES'' must be' ...
        ' numerical date in YYYYMMDD format.'])
    else
        DATEstart = DATES; DATEend = DATES;
    end
elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
        length(num2str(DATES(2)))~=8)
    error(['The value(s) in the second input ''DATES'' must be' ...
        ' numerical date(s) in YYYYMMDD format.'])
else
    DATEstart = DATES(1); DATEend = DATES(2);
end


% Use datenum to accommodate leap years
for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):...
        datenum(num2str(DATEend),'yyyymmdd')

    % Convert back to numerical date
    DATE = str2double(datestr(DATEi','yyyymmdd'));
    % Get default and site specific parameters
    C = getconfig(site,DATE);

    %% Check data and data fields, then load time, range, and snr
    data_stare = loadHALO4bkgCorr(site,DATE,'stare');
    data_ppi   = loadHALO4bkgCorr(site,DATE,'ppi');
    
    %% Sort based on the time stamps
    t_all = [data_stare.time(:); data_ppi.time(:)];
    s_all = [data_stare.snr; data_ppi.snr];
    [~,isort] = sort(t_all);
    snr = s_all(isort,:);
    time_hrs = t_all(isort);
    range_m = data_stare.range(:);
    
    %% Background correction
    % Ripple removal
    snr_corr_1 = correctHaloRipples_v20171018(site,DATE,snr,time_hrs);
    % Correct shape
    [snr_corr_2, ~, ~, ~, ~] = correctBackground_v20171018(...
        snr_corr_1,snr,range_m,time_hrs,'correct_remnant','correct',...
        'ignore',60);
           
    %% Reorganize into original form and write out
    %%--- Stare ---%%
    % Get list of files
    Stare_halo_files = getHALOfileList(site,DATE,'original','stare');
    thedate = num2str(DATE);
    for i = 1:length(Stare_halo_files)
        [Stare_data1,Stare_att1,Stare_dim1] = ...
            load_nc_struct_silent([C.dir_original_stare  ...
            thedate(1:4) '/' Stare_halo_files{i}]);
        
        % Get field names that changed
        corrtd_fnames = checkHALOdata4bkgCorr(site,DATE,'stare');
        
        % Check format and convert into decimal hours
        switch C.time_format
            case 'julian'
                [~,~,~,hh,mm,sss] = jd2date(...
                    Stare_data1.(corrtd_fnames{1})(:));
                original_timehrs = hh(:)+mm(:)/60+sss(:)/3600;
            case 'hours'
                original_timehrs = Stare_data1.(corrtd_fnames{1})(:);
            case 'seconds'
                original_timehrs = Stare_data1.(corrtd_fnames{1})(:)/3600;
        end

        %%--- SIGNAL-TO-NOISE RATIO ---%
        % Assign new values
        % input corrected data according to original orientation
        % Save old data
        Stare_data1.([corrtd_fnames{3} '_uncorrected']) = ...
            Stare_data1.(corrtd_fnames{3});
        % Check dimensions
        [lrow,lcol] = size(Stare_data1.(corrtd_fnames{3}));
        if lrow == length(Stare_data1.(corrtd_fnames{1}))
            Stare_data1.(corrtd_fnames{3}) = snr_corr_2(...
                ismember(time_hrs,original_timehrs),:);
        elseif lcol == length(Stare_data1.(corrtd_fnames{1}))
            Stare_data1.(corrtd_fnames{3}) = snr_corr_2(...
                ismember(time_hrs,original_timehrs),:)';
        else
            error(sprintf(['Length of the field ''time''' ...
                ' doesn''t match either dimensions of the' ...
                ' field ''%s'' in\n%s'], corrtd_fnames{3}, ...
                [C.dir_original_stare thedate(1:4) '/' ...
                Stare_halo_files{i}]))
        end

        % Create attributes
        % Save old data
        Stare_att1.([corrtd_fnames{3} '_uncorrected']) = ...
            Stare_att1.(corrtd_fnames{3});
        if  all(snr == snr_corr_1)
            % bkg*.txt files missing
            Stare_att1.Stare_att1.(corrtd_fnames{3}).long_name = ...
                ['Signal to noise ratio + 1. Background corrected' ...
                ' using Manninen et al. (2016) method.)'];
        else % corrected also by using bkg*.txt files
            Stare_att1.Stare_att1.(corrtd_fnames{3}).long_name = ...
                ['Signal to noise ratio + 1. Background corrected' ...
                ' using Manninen et al. (2016) and Vakkari et al.' ...
                ' (201X) methods.)'];
        end
        if strcmp(site(1:4),'arm-') && ...
                isfield(Stare_att1.global,'datastream')
            %%-- ARM specific --%%
            % add +1 to data level
            dstream0 = Stare_att1.global.datastream;
            d_lev = str2double(Stare_halo_files{i}(length(dstream0)))+1;
            dstream1 = [dstream0(1:strfind(dstream0,'.')+1) ...
                num2str(d_lev)];
            Stare_att1.global.datastream = dstream1;
            
            % generate file name
            fname1 = Stare_halo_files{i};
            fname1(1:length(dstream1)) = dstream1;
        else 
            fname1 = [Stare_halo_files{i}(1:end-length(...
                C.file_format_original_stare)) '_' ...
                C.corrected_filename_identifier ...
                C.file_format_original_stare];
        end
        % write new file
        write_nc_silent([C.dir_corrected_stare thedate(1:4) '/' ...
            fname1], Stare_dim1, Stare_data1, Stare_att1)
    end
    
    
    %%--- PPI ---%%
    % Get list of files
    PPI_halo_files = getHALOfileList(site,DATE,'original','ppi');
    
    for i = 1:length(PPI_halo_files)
        [PPI_data1,PPI_att1,PPI_dim1] = ...
            load_nc_struct_silent([C.dir_original_ppi  ...
            thedate(1:4) '/' PPI_halo_files{i}]);
        
        % Get field names that changed
        corrtd_fnames = checkHALOdata4bkgCorr(site,DATE,'ppi');
        
        % Check format and convert into decimal hours
        switch C.time_format
            case 'julian'
                [~,~,~,hh,mm,sss] = jd2date(...
                    PPI_data1.(corrtd_fnames{1})(:));
                original_timehrs = hh(:)+mm(:)/60+sss(:)/3600;
            case 'hours'
                original_timehrs = PPI_data1.(corrtd_fnames{1})(:);
            case 'seconds'
                original_timehrs = PPI_data1.(corrtd_fnames{1})(:)/3600;
        end

        %%--- SIGNAL-TO-NOISE RATIO ---%
        % Assign new values
        % input corrected data according to original orientation
        % Check dimensions
        % Save old data
        PPI_data1.([corrtd_fnames{3} '_uncorrected']) = ...
            PPI_data1.(corrtd_fnames{3});
        [lrow,lcol] = size(PPI_data1.(corrtd_fnames{3}));
        if lrow == length(PPI_data1.(corrtd_fnames{1}))
            PPI_data1.(corrtd_fnames{3}) = snr_corr_2(...
                ismember(time_hrs,original_timehrs),:);
        elseif lcol == length(PPI_data1.(corrtd_fnames{1}))
            PPI_data1.(corrtd_fnames{3}) = snr_corr_2(...
                ismember(time_hrs,original_timehrs),:)';
        else
            error(sprintf(['Length of the field ''time''' ...
                ' doesn''t match either dimensions of the' ...
                ' field ''%s'' in\n%s'], corrtd_fnames{3}, ...
                [C.dir_original_ppi thedate(1:4) '/' ...
                PPI_halo_files{i}]))
        end

        % Create attributes
        % Save old data
        PPI_att1.([corrtd_fnames{3} '_uncorrected']) = ...
            PPI_att1.(corrtd_fnames{3});
        if  all(snr == snr_corr_1)
            % bkg*.txt files missing
            PPI_att1.(corrtd_fnames{3}).long_name = ...
                ['Signal to noise ratio + 1. Background corrected' ...
                ' using Manninen et al. (2016) method.)'];
        else % corrected also by using bkg*.txt files
            PPI_att1.(corrtd_fnames{3}).long_name = ...
                ['Signal to noise ratio + 1. Background corrected' ...
                ' using Manninen et al. (2016) and Vakkari et al.' ...
                ' (201X) methods.)'];
        end
        if strcmp(site(1:4),'arm-') && ...
                isfield(PPI_att1.global,'datastream')
            %%-- ARM specific --%%
            % add +1 to data level
            dstream0 = PPI_att1.global.datastream;
            d_lev = str2double(PPI_halo_files{i}(length(dstream0)))+1;
            dstream1 = [dstream0(1:strfind(dstream0,'.')+1) ...
                num2str(d_lev)];
            PPI_att1.global.datastream = dstream1;
            
            % generate file name
            fname1 = PPI_halo_files{i};
            fname1(1:length(dstream1)) = dstream1;
        else 
            fname1 = [PPI_halo_files{i}(1:end-length(...
                C.file_format_original_ppi)) '_' ...
                C.corrected_filename_identifier ...
                C.file_format_original_ppi];
        end
        % write new file
        write_nc_silent([C.dir_corrected_ppi thedate(1:4) '/' ...
            fname1], PPI_dim1, PPI_data1, PPI_att1)
    end
end

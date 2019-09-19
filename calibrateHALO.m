function calibrateHALO(site,DATES)
%calibrateHALO loads Halo data, corrects the background, recalculates the
%radial velocity uncertainties, recalculates attenuated backscatter after
%applying focus correction (TBD), and writes data into new files converted
%into Cloudnet naming scheme.
%
% Usage:
% calibrateHALO(site,DATES)
%
% Inputs:
% - site         string, name of the site, e.g. site = 'sodankyla'
% - DATES        scalar or vector, numerical date in decimal hours,
%                e.g. DATES = 20170101 or DATES = [20170101,20170131]
%
% Created 2017-12-11
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
for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):datenum(num2str(DATEend),'yyyymmdd')
    
    % Convert back to numerical date
    DATE = str2double(datestr(DATEi,'yyyymmdd'));
    thedate = datestr(DATEi,'yyyymmdd');
    
    % Get default and site specific parameters
    [C,fnames] = getconfig(site,DATE);
   

    % initialize
    t_all_cell = cell(size(fnames,1),1);
    time_orig_measmode = cell(size(fnames,1),1);
    snr_all_cell = cell(size(fnames,1),1);
    atm_signal_mask_cell = cell(size(fnames,1),1);
    for i_in = 1:size(fnames,1)
        
        % Check path to write to exists
        status = checkHALOpath(site,DATE,'calibrated',fnames{i_in,1},fnames{i_in,2});
        if ~status
            continue
        end
        
        % Get lists of halo files and directory structures
        abc = [fnames{i_in,1} '_' fnames{i_in,2}];
        [loadinfo.(abc).path_to, loadinfo.(abc).files] = getHALOfileList(site,DATE,'original',fnames{i_in,1},fnames{i_in,2});
        if isempty(loadinfo.(abc).files)
            continue
        else
            % Check data and data fields, then load time, range, and snr
            data_tmp = loadHALO4bkgCorr(site,DATE,fnames{i_in,1},fnames{i_in,2});
            if isempty(data_tmp), continue; end
  
            % Check order of time stamps and reorder if needed
            data_tmp = checkHALOtimeStamps(data_tmp);
            
            if isempty(data_tmp)
                continue
            else
                t_all_cell{i_in} = data_tmp.time(:);
                snr_all_cell{i_in} = data_tmp.snr;
                time_orig_measmode{i_in} = zeros(length(data_tmp.time(:)),1)+i_in;
                flag_too_few_profiles = false;
                if size(data_tmp.snr,1)<2
                    flag_too_few_profiles = true;
                    continue
                end
                % Signal filtering, separately for each measurement
                fprintf(['\nFiltering atmospheric signal from %s-%s measurements'...
                    ' for HALO background correction.\nThis might take a while.\n'],fnames{i_in,1},fnames{i_in,2})
                tmp_mask = atmHALOsignalMasking(data_tmp.snr,data_tmp.range,[33 1]);
                atm_signal_mask_cell{i_in} = double(tmp_mask);
            end
        end
    end
    if isempty(cell2mat(snr_all_cell)) || flag_too_few_profiles
        continue
    else
        % convert and sort based on the time stamps
        t_all_mat = cell2mat(t_all_cell);
        time_orig_measmode = cell2mat(time_orig_measmode);
        snr_all_mat = cell2mat(snr_all_cell);
        atm_signal_mask_mat = cell2mat(atm_signal_mask_cell);
        [~,isort] = sort(t_all_mat);
        snr = snr_all_mat(isort,:);
        atm_mask = atm_signal_mask_mat(isort,:);
        time_hrs = t_all_mat(isort);
        range_m = data_tmp.range(:);
        time_orig_measmode = time_orig_measmode(isort);
        
        %% Ripple removal & background correction
        snr_corr_1 = correctHALOripples(site,DATE,snr,time_hrs,range_m);
        if isfield(C,'opt_4_bkg_remnant_profiles') & isfield(C,'num_gates_to_ignore_4bkg_corr')
            snr_corr_2 = correctBackground(snr_corr_1,snr,range_m,time_hrs,...
                'correct_remnant',C.opt_4_bkg_remnant_profiles,...
                'ignore',C.num_gates_to_ignore_4bkg_corr,...
                'cloud_mask',logical(atm_mask));
	else
	    snr_corr_2 = correctBackground(snr_corr_1,snr,range_m,time_hrs,...
	        'cloud_mask',logical(atm_mask));
        end

        
        %% Reorganize into original form and write out
        for i_out = 1:size(fnames,1)
            % Measurement name and type
            abc = [fnames{i_out,1} '_' fnames{i_out,2}];
            
            % Check path to write to exists
            status = checkHALOpath(site,DATE,'calibrated',fnames{i_out,1},fnames{i_out,2});

            if ~status, continue; end
            
            % Get file list
            dir_to_folder_out = getHALOfileList(site,DATE,'calibrated',fnames{i_out,1},fnames{i_out,2});
            
            % file name
            mname = [fnames{i_out,1} '-' fnames{i_out,2}];
            
            % If stare, create one file per day, otherwise keep the original number of files per day
            switch abc
                case {'stare_co','stare_cross'}
                    
                    % load first
                    if isempty(loadinfo.(abc).files), continue; end
                    data0 = load_nc_struct(fullfile(loadinfo.(abc).path_to,loadinfo.(abc).files{1}));
                     
                    % Check order of time stamps and reorder if needed
                    % Check order of time stamps and reorder if needed
                    % Sometimes the last time stamp(s) is(are) in the next day
                    if data0.time(end)<data0.time(end-1)
                        data0.time(find(diff(data0.time)<0+1):end) = data0.time(find(diff(data0.time)<0+1):end) + 24;
                    end                           
 %                   data0 = checkHALOtimeStamps(data0);
                    
                    % Check number of range gates
                    if C.num_range_gates ~= length(data0.(C.field_name_original_range)(:))
                        warning('\n Number of range gates is %s and not %s as specified in the config file --> skipping...',...
                            num2str(length(data0.(C.field_name_original_range)(:))),num2str(C.num_range_gates))
                    else
                        % Convert time into decimal hrs
                        switch C.time_format_original
                            case 'julian'
                                [~,~,~,hh,mm,sss] = jd2date(data0.(C.field_name_original_time)(:));
                                original_time_in_hrs = hh(:)+mm(:)/60+sss(:)/3600;
                            case 'hours'
                                original_time_in_hrs = data0.(C.field_name_original_time)(:);
                            case 'seconds'
                                original_time_in_hrs = data0.(C.field_name_original_time)(:)/3600;
                        end
                        
                        % place corrected signal into the original data struct
                        fname = [C.field_name_original_snr '_' C.corrected_field_name_identifier];
                        [lrow,lcol] = size(data0.(C.field_name_original_snr));
                        if lrow == length(data0.(C.field_name_original_time))
                            data0.(fname) = snr_corr_2(time_orig_measmode==i_out & ismember(time_hrs,original_time_in_hrs),:);
                        elseif lcol == length(data0.(C.field_name_original_time))
                            data0.(fname) = transpose(snr_corr_2(time_orig_measmode==i_out & ismember(time_hrs,original_time_in_hrs),:));
                        end
                        
                        % Convert into cloudnet naming scheme, implement the corrected signal
                        [data1,att1,~] = convert2Cloudnet(site,DATE,fnames{i_out,1},fnames{i_out,2},data0);
                        
                        % If more than one file per day
                        if length(loadinfo.(abc).files)>1
                            
                            % load the rest
                            for i = 2:length(loadinfo.(abc).files)
                                data0 = load_nc_struct(fullfile(loadinfo.(abc).path_to,loadinfo.(abc).files{i}));
                                if C.num_range_gates ~= length(data0.(C.field_name_original_range)(:))
                                    warning('\n Number of range gates is %s and not %s as specified in the config file --> skipping...',...
                                        num2str(length(data0.(C.field_name_original_range)(:))),num2str(C.num_range_gates))
                                else
                                    
                                    % Convert time into decimal hrs
                                    switch C.time_format_original
                                        case 'julian'
                                            [~,~,~,hh,mm,sss] = jd2date(data0.(C.field_name_original_time)(:));
                                            original_time_in_hrs = hh(:)+mm(:)/60+sss(:)/3600;
                                        case 'hours'
                                            original_time_in_hrs = data0.(C.field_name_original_time)(:);
                                        case 'seconds'
                                            original_time_in_hrs = data0.(C.field_name_original_time)(:)/3600;
                                    end
                                    
                                    % place corrected signal into the original data struct
                                    fname = [C.field_name_original_snr '_' C.corrected_field_name_identifier];
                                    [lrow,lcol] = size(data0.(C.field_name_original_snr));
                                    if lrow == length(data0.(C.field_name_original_time))
                                        data0.(fname) = snr_corr_2(time_orig_measmode==i_out & ismember(time_hrs,original_time_in_hrs),:);
                                    elseif lcol == length(data0.(C.field_name_original_time))
                                        data0.(fname) = transpose(snr_corr_2(time_orig_measmode==i_out & ismember(time_hrs,original_time_in_hrs),:));
                                    end
                                                                
                                    % Convert into cloudnet naming scheme, implement the corrected signal
                                    tmpdata1 = convert2Cloudnet(site,DATE,fnames{i_out,1},fnames{i_out,2},data0);
                                    
                                    % Stack on the first files fields that are not scalars AND not 'range'
                                    fnamesdata1 = fieldnames(tmpdata1);
                                    for i2 = 1:length(fnamesdata1)
                                        % Stack only those fields where time is one dimension
                                        if ~any(size(tmpdata1.(fnamesdata1{i2}))==size(tmpdata1.time(:),1)), continue; end
                                        % Don''t stack range
                                        if strcmp(fnamesdata1{i2},'range'), continue; end
                                        % Stack the rest...
                                        data1.(fnamesdata1{i2}) = [data1.(fnamesdata1{i2}); tmpdata1.(fnamesdata1{i2})];
                                    end
                                end
                            end
                        end
                        
                        % Sometimes the last time stamp is from the next day
                        if data1.time(end)<data1.time(end-1)
                            data1.time(end) = data1.time(end)+24;
                        end
                        
                        dim1 = struct('time',length(data1.time),'range',length(data1.range));
                        fndate = loadinfo.(abc).files{1}(1:8);
                        if ~strcmp(fndate,thedate)
			    fndate = thedate;
                        end
                        
                        % correct focus TBD
                        data1 = correctHALOfocus(site,DATE,abc,data1);
                        
                        % write new file
                        write_nc_struct(fullfile([dir_to_folder_out '/' fndate '_' site '_halo-doppler-lidar-' num2str(C.halo_unit_id) '-' mname '.nc']),dim1,data1,att1);
                    end
                otherwise
                    for i = 1:length(loadinfo.(abc).files)
                        %-- load --%
                        if isempty(loadinfo.(abc).files), continue; end
                        data0 = load_nc_struct(fullfile(loadinfo.(abc).path_to,loadinfo.(abc).files{i}));

                        % Check order of time stamps and reorder if needed
                        % Sometimes the last time stamp(s) is(are) in the next day
                        if data0.time(end)<data0.time(end-1)
			    data0.time(find(diff(data0.time)<0+1):end) = data0.time(find(diff(data0.time)<0+1):end) + 24;
                        end                           
%                        data0 = checkHALOtimeStamps(data0);
                        
                        % Check number of range gates
                        if C.num_range_gates ~= length(data0.(C.field_name_original_range)(:))
                            warning('\n Number of range gates is %s and not %s as specified in the config file --> skipping...',...
                                num2str(length(data0.(C.field_name_original_range)(:))),num2str(C.num_range_gates))
                        else
                            
                            % Convert time into decimal hrs
                            switch C.time_format_original
                                case 'julian'
                                    [~,~,~,hh,mm,sss] = jd2date(data0.(C.field_name_original_time)(:));
                                    original_time_in_hrs = hh(:)+mm(:)/60+sss(:)/3600;
                                case 'hours'
                                    original_time_in_hrs = data0.(C.field_name_original_time)(:);
                                case 'seconds'
                                    original_time_in_hrs = data0.(C.field_name_original_time)(:)/3600;
                            end
                            
                            % place corrected signal into the original data struct
                            fname = [C.field_name_original_snr '_' C.corrected_field_name_identifier];
                            [lrow,lcol] = size(data0.(C.field_name_original_snr));
                            if lrow == length(data0.(C.field_name_original_time))
                                if size(data0.(C.field_name_original_snr),1)<2, continue; end
                                data0.(fname) = snr_corr_2(time_orig_measmode==i_out & ismember(time_hrs,original_time_in_hrs),:);
                            elseif lcol == length(data0.(C.field_name_original_time))
                                if size(data0.(C.field_name_original_snr),2)<2, continue; end
                                data0.(fname) = transpose(snr_corr_2(time_orig_measmode==i_out & ismember(time_hrs,original_time_in_hrs),:));
                            end
                            
                            % Convert into cloudnet naming scheme, implement the corrected signal
                            [data1,att1] = convert2Cloudnet(site,DATE,fnames{i_out,1},fnames{i_out,2},data0);
                            
                            % If more than one file per day
                            check_time = [thedate '_' datestr(decimal2daten(data0.time(1),datenum(thedate,'yyyymmdd')),'HHMMSS')];
                            if length(loadinfo.(abc).files)>1                
                                fndate = loadinfo.(abc).files{i}(1:15);
                                if ~strcmp(fndate,check_time)
			            fndate = check_time;
                                end
                            else
                                fndate = loadinfo.(abc).files{i}(1:8);
                                if ~strcmp(fndate,thedate)
   			            fndate = thedate;
                                end
                            end
                            
                            %% Sometimes the last time stamp is from the next day
                            %if data1.time(end)<data1.time(end-1)
                            %    data1.time(end) = data1.time(end)+24;
                            %end
                            
                            % correct focus
                            data1 = correctHALOfocus(site,DATE,abc,data1);
                            
                            % Dimensions
                            dim1 = struct('time',length(data1.time),'range',length(data1.range));
                            
                            % write new file
                            write_nc_struct(fullfile([dir_to_folder_out '/' fndate '_' site '_halo-doppler-lidar-' num2str(C.halo_unit_id) '-' mname '.nc']),dim1,data1,att1);
                        end
                    end
            end
        end
    end
end
end

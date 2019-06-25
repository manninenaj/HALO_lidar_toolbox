function [snr1,steps] = correctHALOripples(site,DATE,snr0,t_snr)   
%correctHALOripples corrects the ripples in the HALO Doppler lidar
% background, see Vakkari et al. (2019).
%
% Usage:
% snr1 = correctHALOripples(site,DATE,snr0,t_snr)
% [snr1,steps] = correctHALOripples(site,DATE,snr0,t_snr)
%
% Inputs:
% - site     string, name of the site, e.g. 'kuopio'
% - DATE     scalar, numerical date, e.g. 20171231
% - snr0     array (time,height), uncalibrated SNR profiles
% - t_snr    array (time,1), timestamps for SNR profiles,
%
% Outputs:
% - snr1     array (time,height), ripple corrected SNR
% - steps    array, locations of expected steps in SNR
%
% Modified 2019-06-25
% Antti Manninen
% antti.manninen(at)fmi.fi
% Finnish Meteorological Institute

C = getconfig(site,DATE);
% In case bkg directory is not specified, skip ripple correction
if ~isfield(C,'dir_background_txt') || ...
        length(C.dir_background_txt)<=2
    snr1 = snr0;
    steps = [];
    return
end

% Set dates
daten = datenum(num2str(DATE),'yyyymmdd');
time_snr_dnum = decimal2daten(t_snr,daten);
thedate = num2str(DATE);

% Get path and list of files
[bkg_path, files_bkg] = getHALOfileList(site, DATE, 'background', ...
                                            'txt');
bkg_path_tday = bkg_path; % for later use
files_bkg_tday = files_bkg; % for later use

bkg_time1 = nan; % initialize
if ~isempty(files_bkg)
    % Get the timestamp of the first bkg file from today
    bkg_time1 = datenum([files_bkg{1}(12:17) files_bkg{1}(19:24)], ...
                        'ddmmyyHHMMSS');
end

% Go back in time to get the timestamp of the last bkg profile
% Hopefully the last bkg file is not too far back in time...
fprintf("Looking for the last bkg profile recorded before '%s'.\n", ...
        datestr(time_snr_dnum(1),'yyyy-mm-dd HH:MM:SS'));
bkg_time_last = bkg_time1;
files_bkg_yday = [];
daten_yday = daten;
while isnan(bkg_time1) | bkg_time_last > time_snr_dnum(1)
    daten_yday = daten_yday-1;
    DATE_yday = str2num(datestr(daten_yday,'yyyymmdd'));
    if DATE_yday >= C.parameters_valid_from_including
        fprintf("Checking '%s'.\n", datestr(daten_yday,'yyyy-mm-dd HH:MM:SS'));
        [bkg_path_yday, files_bkg_yday] = getHALOfileList(site, ...
                                                          DATE_yday, ...
                                                          'background', ...
                                                          'txt');
        if isempty(files_bkg_yday)
            continue
        else
            % Check the last bkg profile
            bkg_time_last = datenum([files_bkg{end}(12:17) ...
                                files_bkg{end}(19:24)], 'ddmmyyHHMMSS');
        end
    else
        if isempty(files_bkg)
            snr1 = snr0;
            steps = [];
            msg_1 = ["No background files found for site '%s' that" ...
                     " exist before '%s'.\n"];
            fprintf(msg_1,site,thedate);
            return
        else
            msg_2 = "No older bkg files found than '%s'.\n";
            fprintf(msg_2,[bkg_path files_bkg{1}]);
            break
        end
    end
end

% Skip ripple removal if no bkg files for the site
if isempty(files_bkg_tday) && isempty(files_bkg_yday)
    snr1 = snr0;
    steps = [];
    return
end

% Calculate P_bkg and P_fit for todays and/or previous bkg profiles
% Initialize
P_bkg_yday = [];
P_fit_yday = [];
bkg_times_yday = [];
P_bkg_tday = [];
P_fit_tday = [];
bkg_times_tday = [];
if ~isempty(files_bkg_yday)
    [P_bkg_yday, P_fit_yday, bkg_times_yday] = calculateBKGtxt(bkg_path_yday, ...
                                                files_bkg_yday, daten, ...
                                                C.num_range_gates);    
end
if ~isemtpy(files_bkg_tday)
    [P_bkg_tday, P_fit_tday, bkg_times_tday] = calculateBKGtxt(bkg_path_tday, ...
                                                files_bkg_tday, daten, ...
                                                C.num_range_gates);
end

% Put all together
P_bkg = [P_bkg_yday(:),P_bkg_tday(:)];
P_fit = [P_fit_yday(:),P_fit_tday(:)];
bkg_times = [bkg_times_yday(:);bkg_times_tday(:)];

% Scale to the snr time resolution, assing matching bkg profile
P_bkg_scl = nan(size(snr0));
P_fit_scl = nan(size(snr0));
for ib = 1:length(bkg_times)
    P_bkg_scl(find(bkg_times(ib)<=time_snr_dnum),:) = P_kbg(ib,:);
    P_fit_scl(find(bkg_times(ib)<=time_snr_dnum),:) = P_fit(ib,:);
end

% Remove bkg profiles that contain only nans
P_bkg(all(isnan(P_bkg),2),:) = 0; 
P_fit(all(isnan(P_fit),2),:) = 0;
bkg_times(all(isnan(P_fit),2)) = nan;

% Find steps
steps = nan(1,length(bkg_times)-1);
for i = 2:length(bkg_times)
    if ~isempty(find(time_snr_dnum>bkg_times(i),1,'first'))
        steps(i-1) = find(time_snr_dnum>bkg_times(i),1, ...
                                   'first');
    end    
end

% For now, P_amp is calculated only for few test sites
if isfield(C,'dir_housekeeping')
    path_to_P_amp = [C.dir_housekeeping thedate(1:4) ...
                     '_amp_ave_resp.mat'];
    if exist(path_to_P_amp)==2
        P_amp = load(path_to_P_amp);
        P_amp = P_amp.P_amp;
       P_amp_exists = true;
    else
        P_amp_exists = false;
    end
else
    P_amp_exists = false;
end

% Remove ripples by using P_fit (and P_amp if available)
istep = 1;
snr1 = nan(size(snr0));
for i = 1:size(snr0,1)
    if istep < size(P_bkg,1)
        if i >= steps(istep)
            istep = istep + 1;
        end        
    end
    if P_amp_exists
        P_noise = P_fit + repmat(transpose(P_amp),size(P_fit,1),1);
    else
        P_noise = P_fit;
    end
    % In case some SNR profiles cannot be associated with bkg
    % profiles, set P_bkg and P_fit to '1' so SNR will not be
    % affected in those cases 
    P_noise(P_noise==0)= 1;
    P_bkg(P_bkg==0) = 1;
    % Remove ripples
    snr1(i,:) = snr0(i,:) .* (P_bkg(istep,:) ./ P_noise(istep,:));
    end
end
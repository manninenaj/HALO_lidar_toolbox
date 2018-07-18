function [snr1,step_locations] = correctHALOripples(site,DATE,...
    snr0,t_snr)
%correctHALOripples corrects the ripples in the HALO bakcground, see
%Vakkari et al. (201?)

C = getconfig(site,DATE);
daten = datenum(num2str(DATE),'yyyymmdd');

% Get the date folder structure
[bkg_path, files_bkg] = getHALOfileList(site,DATE,'background','txt');

if ~isempty(files_bkg)
% if exist(path_bkg,'dir') == 7 % if bkg files exist
    [b_file, b_fit, bkg_times] = calculateBKGtxt(bkg_path,daten,...
        size(snr0,2));
%     b_file(:,1:3) = nan; b_fit(:,1:3) = nan;
    b_file(all(isnan(b_file),2),:) = []; 
    b_fit(all(isnan(b_fit),2),:) = [];
    
    % Find steps
    [time_snr_dnum] = decimal2daten(t_snr,daten);
    step_locations = nan(1,length(bkg_times)-1);
    for i = 2:length(bkg_times)
        if ~isempty(find(time_snr_dnum>bkg_times(i),1,'first'))
            step_locations(i-1) = find(time_snr_dnum>bkg_times(i),1,...
                'first');
        end
    end
    
    % Correct ripples
    istep = 1;
    snr1 = nan(size(snr0));
    for i = 1:size(snr0,1)
        if istep < size(b_file,1)
            if i >= step_locations(istep)
                istep = istep + 1;
            end
        end
        b_snr = b_file(istep,:)./b_fit(istep,:);
        snr1(i,:) = snr0(i,:) + b_snr - 1;
    end
else 
    snr1 = snr0;
    step_locations = [];
end
end


function estimateHALOaverageAmplifierResponse(site,YEAR)
% Still under development...

if ~strcmp(site,'juelich')
    return
else
    P_bkg_res = cell(1,1);
    i = 1;
    for daten = datenum(YEAR,01,01):datenum(YEAR,12,31)
        thedate = datestr(daten,'yyyymmdd');
        DATE = str2num(thedate);
        C = getconfig(site,DATE);
        % Get path
        bkg_path = getHALOfileList(site,DATE,'background','txt');
        % Load the background profiles
        [P_bkg, P_fit] = calculateBKGtxt(bkg_path,daten,C.num_range_gates);
        P_bkg_res{i} = P_bkg - P_fit;
        i = i + 1;
    end   
    P_bkg_res_all = cell2mat(P_bkg_res(:));
    if license('test','wavelet_toolbox')
        P_amp = wdenoise(nanmean(P_bkg_res_all),'Wavelet','sym8');
        P_amp = P_amp(:);
        save([C.dir_housekeeping num2str(YEAR) '_amp_ave_resp.mat'],'P_amp')
    else
        warning('Wavelet toolbox missing, P_amp estimation skipped.')
    end
end
end


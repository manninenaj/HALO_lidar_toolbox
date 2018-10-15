function calculateHALObkgMedianWaveStructure(site)
%calculateHALObkgWaveStructure calculates the median values for the 
%persistent structure in the HALO background, which is caused by the 
%amplifier's response to the transmitted pulse (Vakkari et al. 2018).
%
% Usage:
% calculateHALObkgWaveStructure(site)
%
% Inputs:
% - site         string, name of the site, e.g. site = 'sodankyla'
%
% Created 2018-10-02
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% INAR
% University of Helsinki, Finland

% Check inputs
if ~ischar(site)
    error('The first input ''site'' must be a string.')
end

% These should go into the database!!
DATES = {'20180124_010000-20180126_120000';...
         '20180201_180000-20180205_090000';...
         '20180213_000000-20180214_140000';...
         '20180315_030000-20180323_120000'};

     
     
for iYear = 2010:2090

    thedate = [num2str(iYear) '0125'];
    DATE = str2num(thedate);
    C = getconfig(site,DATE);
    
    bkg_path = C.dir_background_txt;
    bkg_path = strrep(bkg_path,'+YYYY+',thedate(1:4));
    bkg_path = strrep(bkg_path,'+MM+',thedate(5:6));
    bkg_path = strrep(bkg_path,'+DD+',thedate(7:8));
    
    [dirto,files] = getHALOfileList(site,DATE,'original','stare','co');
    [~,~,dims] = load_nc_struct(fullfile([dirto files{1}]),{'time','range'});
    [bkg_out, fit_out, bkg_times] = calculateBKGtxt(bkg_path,datenum(thedate,'yyyymmdd'),dims.range)
        
        
end

end

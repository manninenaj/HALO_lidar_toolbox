function aftercorrectHALObkgWaves(site,DATES)
%aftercorrectHALObkgWaves corrects for the wave-like forms in the signal
%after the calibration, i.e. bkg correction and focus correction. It reads
%predefined period of data where these wave-like shapes can be extracted
%without contamination of atmospheric signal. If such a period is not
%defined for the site and when the same HALO instrument was deployed with
%the same measurement parameters, the functions returns a warning and no
%wave correction is carried out. Then, signal is corrected, beta values
%re-calculated and a new netcdf file is written.
%
% Usage:
% aftercorrectHALObkgWaves(site,DATES)
%
% Inputs:
% - site         string, name of the site, e.g. site = 'sodankyla'
% - DATES        scalar or vector, numerical date in decimal hours,
%                e.g. DATES = 20170101 or DATES = [20170101,20170131]
%
% Created 2018-10-02
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% INAR
% University of Helsinki, Finland

 Check inputs
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
end

end


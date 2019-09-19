function calculateHALOsigma2vadProduct(site,DATES,elevangle)
%calculateHALOsigma2vadProduct reads ppi files and calculates... 
%
% Usage:
%  calculateHALOsigma2vadProduct(site,DATES,elevangle)
%
% Inputs:
% -site        String, site name, e.g. site = 'kuopio'
% -DATES       Scalar or vector, numeric, e.g. DATES = 20170401 
%              or DATES = [20170401 20170431]
% -elevangle   string, elevation angle 0-90
%
% Created 2019-09-03
% Antti Manninen
% Finnish Meteorological Institute
% 

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
elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
        length(num2str(DATES(2)))~=8)
    error(['The value(s) in the second input ''DATES'' must be' ...
        ' numerical date(s) in YYYYMMDD format.'])
else
    DATEstart = DATES(1); DATEend = DATES(2);
end
if ~ischar(elevangle) || length(elevangle) ~= 2 || (~isempty(str2num(elevangle)) && str2num(elevangle)<0 || str2num(elevangle)>90) 
    error('The 3rd input must be a string and no longer than 2 characters specifying the elevation angle 0-90 degrees.')
end

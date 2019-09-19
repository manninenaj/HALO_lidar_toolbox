function yd = my_doy(thedate)
% DOY Ordinal number of day in a year.
%
%  DOY = MY_DOY(YYYYMMDD) returns the ordinal day number in the given year 
%
%  Returns value in the same format (string or double) as input


output = 'string';

if isnumeric(thedate)
  thedate = num2str(thedate);
  output = 'numeric';
end

if length(thedate) ~= 8
  error('Error: Date should be in format YYYYMMDD' );
  return
end

yd = datenum(thedate, 'yyyymmdd') - datenum(thedate(1:4), 'yyyy') + 1;

if strcmp(output,'string')
  yd = num2str(yd);
end



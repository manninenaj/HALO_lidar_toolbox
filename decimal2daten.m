function [t_dnum] = decimal2daten(decimal_time,date_dnum)
%DECIMAL2DATEN converts decimal time to matlab serial date format
%
% Usage: [t_dnum] = decimal2daten(decimal_time,date_dnum)
%
% Inputs:
% - decimal_time      decimal time vector (e.g. 0.01 0.02 ... 23.98 23.99)
% - date_dnum         date in matlab serial date format
%
% Output:
% - t_dnum            time vector in matlab serial date format
% 
% 2016 July 17
% Antti Manninen
% antti.j.manninen@helsinki.fi

if nargin < 2
    error 'Not enough input arguments!'
elseif nargin == 2
    switch isvector(decimal_time) && ...
            sum(~isfinite(decimal_time)) == 0 && ...
            sum(decimal_time < 0) == 0
        case 1
            switch isscalar(date_dnum) && ...
                    date_dnum > 0 && ...
                    isfinite(date_dnum)
                case 1
                    decimal_time = decimal_time(:);
                    decimal_time = double(decimal_time);
                    Hrs = floor(decimal_time);
                    Min = floor(((decimal_time) - Hrs) * 60);
                    Sec = floor(((((decimal_time) - Hrs) * 60) - Min) * ...
                        60);
                    t_dnum = repmat(date_dnum,size(Hrs)) + ...
                        datenum(0,0,0,Hrs,Min,Sec);
                otherwise
                    error 'Second input must be a positive finite scalar!'
            end
        otherwise
            error 'First input must contain only positive finite values!'
    end     
else
    error 'Too many input arguments!'
end

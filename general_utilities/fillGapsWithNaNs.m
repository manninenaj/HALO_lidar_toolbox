function [ X_out,t_out ] = fillGapsWithNaNs(t_in,X_in,dt)
%fillGapsWithNaNs looks for inconsecutive time stamps and adds a nan in between,
%right after the previous time stamp
%
% Usage:
% [X_out,t_out] = fillGapsWithNaNs(t_in,X_in)
% [X_out,t_out] = fillGapsWithNaNs(t_in,X_in,dt)
%
% Inputs:
% - t      Column vector containing the time stamps
% - X_in   Input array in which the nans are to be added
% - dt     By default the function determines the inconsecutiveness by
%           looking at the median difference in time stamps. The 'dt',
%           which is by default 1, adds an adjustment to the allowed
%           difference: abs(t(i)-t(i-1)) > median(diff(t)) * dt 
%
% Outputs:
% - X_out  Output array with nans added in between inconsecutive 
%           time stamps
% - t_out  Output time stamps with nans added in between inconsecutive 
%           time stamps
%
% 2016-01-14
% Antti Manninen
% antti.j.manninen@helsinki.fi

% Check the number of inputs and set the defaults
if nargin < 2
    error 'Not enough inputs! Type > help fillGaps <.'
elseif nargin < 3
    dt = 1;
end

% Check inputs
% validateattributes(t_in,{'double'},{'2d','column','increasing','nonempty'},'','',1)
% validateattributes(X_in,{'double'},{'2d','nrows',length(t_in)},'','',2)
% validateattributes(dt,{'double'},{'scalar','positive'},'','',3)

X_out(1,:) = X_in(1,:);
t_out(1) = t_in(1);
i_out = 2;
i_in = 2;
while i_in <= size(X_in,1)
    switch abs(t_in(i_in) - t_in(i_in - 1)) > nanmedian(diff(t_in)) * dt
        case 0
            X_out(i_out,:) = X_in(i_in,:);
            t_out(i_out) = t_in(i_in);
            i_in = i_in + 1;
            i_out = i_out + 1;
        case 1
            t_out(i_out) = t_in(i_in-1)+nanmedian(diff(t_in)) * dt;
            X_out(i_out,:) = nan(1,size(X_in,2));
            i_out = i_out + 1;
            X_out(i_out,:) = X_in(i_in,:);
            t_out(i_out) = t_in(i_in);
            i_in = i_in + 1;
            i_out = i_out + 1;
    end
end
t_out = t_out(:);
end


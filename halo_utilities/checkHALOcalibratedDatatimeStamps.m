function [data] = checkHALOcalibratedDatatimeStamps(data)
%checkHALOcalibratedDatatimeStamps sort time stamps and sorts all the
%related fields according to the sorted order.
%
% Usage: data = checkHALOcalibratedDatatimeStamps(data)
%
% Inputs:
% - data        struct
%
% Outputs:
% - data        struct
%
%
% Created 2018-09-13
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% INAR
% University of Helsinki, Finland

if ~isstruct(data)
    error('Input must be a struct!')
end

% Sort time stamps and get order
[~,iorder] = sort(data.time(:));
% Get fieldnames
fnames = fieldnames(data);
for ifn = 1:length(fnames)
    % Sort those fields in which time is one dimension
    if ~any(size(data.(fnames{ifn}))==size(data.time(:),1)), continue; end
    % Check the orientation, and sort
    [lrow,lcol] = size(data.signal);
    if lrow == length(data.time)
        data.(fnames{ifn}) = data.(fnames{ifn})(iorder,:);
    elseif lcol == length(data.time)
        data.(fnames{ifn}) = data.(fnames{ifn})(:,iorder);
    end
end


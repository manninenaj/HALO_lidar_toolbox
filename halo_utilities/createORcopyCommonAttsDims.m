function [data_out,att_out,dim_out] = createORcopyCommonAttsDims(data_in,C,att_in)
%createORcopyCommonAttsDims creates common attirubutes for all Doppler lidar products
%                 
%         - - - - - - - - - - - - - - - - - - 
%          \     |               |   |
%      range\    |range          |   |
%            \   |               |   |
%             \  |     height_agl|   |height
%              \ |               |   |  
%               ____             |   |
%              _|__|_ - - - - - -|- -|- - -
%        lidar |    |            |   |   |
%           ___|____|___         |   |   |
%           | __ __ __ |         |   |   |
%           | || || || |         |   |   |
%   building| __ __ __ |         |   |   |
%           | || || || |         |   |   |altitude_instrument 
%           | __ __ __ |         |   |   |
%           | || || || |         |   |   |
%      _____|__________|_________|_  |   | 
%      ground level          |       |   |                           
%                    altitude|       |   |  
%                            |       |   |  
%                           ~~~~~~~~~~~~~~~~~~~ 
%                           mean sea level
%
% Created 2019-09-17
% Antti Manninen
% Finnish Meteorological Institute
% antti.manninen@fmi.fi

if ~isfield(data_in,'elevation')
    elev = 90; % Assume vertical profile data
else 
    elev = data_in.elevation;
end
data_out.elevation = elev;

% Add dims
if ~isfield(data_in,'range')
    dim_out.range = length(data_in.height);
else
    dim_out.range = length(data_in.range);
end

% NOTE: With old data, height = range and assumes instrument at ground level
% Add range if it does not exist, possible with products (old version)
% where products were always calculated w.r.t. height agl
if ~isfield(data_in,'range')
    data_in.range = data_in.height;
    data_out.range = data_in.height;
    att_out.range = create_attributes(...
        {'range'},...
        'Range (vertical) from the instrument', ...
        'm',...
        [],...
        'Towards zenith.');
else
    data_out.range = data_in.range;
    att_out.range = create_attributes(...
        {'range'},...
        'Range of measurement from the instrument', ...
        'm',...
        [],...
        'Range towards the line-of-sight.');
end

att_out.range.axis = 'Z';

% heights
if isfield(C,'altitude_instrument_level_m_asl') && isfield(C,'altitude_ground_level_m_asl')
    actual_height_above_ground = C.altitude_instrument_level_m_asl - C.altitude_ground_level_m_asl;
    actual_instrument_altitude_asl = C.altitude_instrument_level_m_asl;
    actual_site_altitude_asl = C.altitude_ground_level_m_asl;
    cmnt = '.';
else
    actual_height_above_ground = 0; % assumes instrument at ground level
    actual_instrument_altitude_asl = C.altitude_in_meters;
    cmnt = ', assumes instrument at ground level becuase no information given.';
    actual_site_altitude_asl = C.altitude_in_meters;
end
%if length(unique(elev))~=1 && any(diff(elev)>.5)
%    for i = 1:length(elev)
%        data_out.height_agl(i,:) = data_in.range .* sind(elev(i)) + actual_height_above_ground;
%        data_out.height_asl(i,:) = data_in.range .* sind(elev(i)) + actual_instrument_altitude_asl;
%    end
% assume the elevation was supposed to be the same...
%else
data_out.height_agl = data_in.range .* sind(nanmedian(elev)) + actual_height_above_ground;
data_out.height_asl = data_in.range .* sind(nanmedian(elev)) + actual_instrument_altitude_asl;
%end
att_out.height_agl = create_attributes(...
    {'range'},...
    'Height above ground level', ...
    'm',...
    [],...
    ['Height_agl = range * sin(elevation) + height of instrument above ground' cmnt]);
att_out.height_agl.axis = 'Z';
att_out.height_asl = create_attributes(...
    {'range'},...
    'Height above mean sea level', ...
    'm',...
    [],...
    ['Height_asl = range * sin(elevation) + altitude of instrument above mean sea level' cmnt]);
att_out.height_asl.axis = 'Z';

% altitude
data_out.altitude_site = actual_site_altitude_asl;
att_out.altitude_site = create_attributes(...
    {},...
    ['Altitude of site above mean sea level' cmnt], ...
    'm');
% altitude of instrument
data_out.altitude_instrument = actual_instrument_altitude_asl;
att_out.altitude_instrument = create_attributes(...
    {},...
    ['Altitude of instrument above mean sea level' cmnt], ...
    'm');

% latitude
if ~isfield(data_in,'latitude')
    if ~isfield(C,'latitude')
        data_out.latitude = 0;
    else
        data_out.latitude = C.latitude;
    end
 else 
    data_out.latitude = data_in.latitude;
end
att_out.latitude = create_attributes(...
    {},...
    'Latitude of lidar', ...
    'degrees_north');
att_out.latitude.standard_name = 'latitude';

% longitude
if ~isfield(data_in,'longitude')
    if ~isfield(C,'longitude')
        data_out.longitude = 0;
    else
        data_out.longitude = C.longitude;
    end
 else 
    data_out.longitude = data_in.longitude;
end
att_out.longitude = create_attributes(...
    {},...
    'Longitude of lidar', ...
    'degrees_east');
att_out.longitude.standard_name = 'longitude';

% Order fields
data_out = orderfields(data_out);
att_out  = orderfields(att_out);
dim_out = orderfields(dim_out);



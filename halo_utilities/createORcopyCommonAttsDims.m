function [data_out,att_out,dim_out] = createORcopyCommonAttsDims(data_in,C)
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
%  building | __ __ __ |         |   |   |
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


% Add dims
dim_out.range = length(data_in.height);

% NOTE: With old data, height = range and assumes instrument at ground level
if ~isfield(data_in,'height_agl') % check from input 'data_in' 
        % If add range if it does not exist
    	if ~isfield(data_in,'range')
            data_out.range = data_in.height;
        end
        att_out.range = create_attributes(...
            {'range'},...
            'Range from instrument', ...
    	    'm');
        att_out.range.axis = 'Z';
        % height agl
	data_out.height_agl = data_in.height;
	att_out.height_agl = create_attributes(...
    	    {'range'},...
	    'Height above ground level', ...
	    'm',...
	    [],...
	    ['Height_agl = range (vertically pointing), assumes instrument is at ground level']);
        att_out.height_agl.axis = 'Z';

        
        % height
	att_out.height = create_attributes(...
	    {'range'},...
	    'Height above mean sea level', ...
	    'm',...
	    [],...
	    ['Height above mean sea level']);
	att_out.height.axis = 'Z';

        % If height_agl did not exist, but if instruments and ground altitudes are specified 
        % in halo_config.txt, correct height_agl
	if isfield(C,'altitude_ground_level_m_asl') && isfield(C,'altitude_instrument_level_m_asl') % check from config
	    % create height_agl
	    actual_height_above_ground = C.altitude_instrument_level_m_asl - C.altitude_ground_level_m_asl;
	    data_out.height_agl = data_in.height + actual_height_above_ground;
	    att_out.height_agl.comment = 'Height above ground level';

	    % old height = range --> correct to height asl
	    data_out.height = data_in.height + C.altitude_instrument_level_m_asl;

            % altitude
            data_out.altitude = C.altitude_ground_level_m_asl;
            att_out.altitude = create_attributes(...
                {},...
                'Altitude of site above mean sea level', ...
                'm');
            % altitude of instrument
            data_out.altitude_instrument = C.altitude_instrument_level_m_asl;
            att_out.altitude_instrument = create_attributes(...
                {},...
                'Altitude of instrument above mean sea level', ...
                'm');

	else         
            % Leave height_agl as is, added for brevity (but commented out)
    	    %% data.height_agl = data.height_agl;
 
	    % Use old altitude in meters parameter for correcting height asl, assume instrument at ground
	    data_out.height = data_in.height + C.altitude_in_meters;
	    att_out.height.comment = 'Height above mean sea level, assumes instrument is at ground level';

            % altitude
            data_out.altitude = C.altitude_in_meters;
            att_out.altitude = create_attributes(...
                {},...
                'Altitude of site above mean sea level', ...
                'm');

            % altitude of instrument, assume at ground level, so the same as altitude
            data_out.altitude_instrument = data_in.altitude;
            att_out.altitude_instrument = att_in.altitude;
	end
end

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



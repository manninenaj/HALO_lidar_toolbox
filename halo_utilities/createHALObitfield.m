function [ bitfield ] = createHALObitfield(th,time,beta,...
    aero_top,flux,epsilon,shear_vector,sun_rise_set,fubarfield,precip_bit)
%createHALObitfield creates a bitfield for boundary layer
%classification
bitfield = nan(size(beta));
for i = 1:size(beta,1)
    for j = 1:size(beta,2)
        
        % Aerosol layer bit
        switch j <= aero_top(i)
            case 0
                bit_layer = 0;
            case 1
                bit_layer = 1;
        end
        
        if ~isempty(flux) % Is flux measurement is available
            % Heat flux bit
            switch flux(i) > th.heatflux
                case 0
                    bit_flux = 0;
                case 1
                    bit_flux = 1;
            end
            % Heat flux neutral
            switch abs(flux(i) - 0) > th.heatflux_hi
                case 0
                    bit_flux_hi = 0;
                case 1
                    bit_flux_hi = 1;
            end
        else
            % Heat flux neutral
            switch time(i) > sun_rise_set(1) && time(i) < sun_rise_set(2)
                case 0
                    bit_flux = 0;
                case 1
                    bit_flux = 1;
            end
            bit_flux_hi = 1; % not flux meas., only stable/unstable classes
        end
        
        % Dissipation rate bit
        switch real(log10(epsilon(i,j))) > th.epsilon
            case 0
                bit_eps = 0;
            case 1
                bit_eps = 1;
        end
        
        %  Wind shear bit
        switch shear_vector(i,j) > th.windshear
            case 0
                bit_shear = 0;
            case 1
                bit_shear = 1;
        end
        
        % In contact with the surface
        switch fubarfield(i,j) == 2
            case 0
                bit_contact_w_s = 0;
            case 1
                bit_contact_w_s = 1;
        end
        
        % Convective
        switch fubarfield(i,j) == 6
            case 0
                bit_convective = 0;
            case 1
                bit_convective = 1;
        end
        
        % In cloud bit
        switch fubarfield(i,j) == 4
            case 0
                bit_incloud = 0;
            case 1
                bit_incloud = 1;
        end
        
        % Cloud driven bit
        switch fubarfield(i,j) == 3
            case 0
                bit_cloud_driven = 0;
            case 1
                bit_cloud_driven = 1;
        end
        
        % precipitation bit
        switch precip_bit(i)
            case 0
                bit_precipitation = 0;
            case 1
                bit_precipitation = 1;
        end
        
        bitfield(i,j) = uint16(bit_layer + ...
            2*bit_eps + ...
            4*bit_contact_w_s + ...
            8*bit_flux + ...
            16*bit_flux_hi + ...
            32*bit_shear + ...
            64*bit_incloud + ...
            128*bit_cloud_driven + ...
            256*bit_convective + ...
            512*bit_precipitation);
    end
end
end

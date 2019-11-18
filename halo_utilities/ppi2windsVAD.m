    function [S_out] = ppi2windsVAD(S_in,varargin)
%ppi2windsVAD calculates wind retrieval and uncertainties based on methods
% described by: Paeschke et al. (2015) and Newsom et al. (2017).
%
% Inputs (units)[dimensions]
% - S_in              struct with the followind fields:
%    .time            vector, time in numerical format
%    .range           vector, range from the instrument
%    .azimuth         vector, azimuth angles (degrees 0-360)[time 1]
%    .elevation       vector, elevation angles (degrees 0-360)[time 1]
%    .velocity_raw    matrix, radial velocities (m s-1)[time range]
%    .velocity        matrix, noise filtered radial velocities (m s-1)[time range]
%    .snr             matrix, signal-to-noise ratio [time range]
%    .velocity_error  matrix, radial velocity uncertainties (m s-1)[time range]
%    varargin:        parameter-value pairs
%      - fit_error     'fit_error',false, true is default
%
% Outputs (units)[dimensions]:
% - S_out            struct with the followind fields:
%    .u              u-wind (m s-1)[range]
%    .v              v-wind (m s-1)[range]
%    .w              w-wind (m s-1)[range]
%    .ws             wind speed (m s-1)[range]
%    .wd             wind direction (degrees)[range]
%    .u_error        u-wind error (m s-1)[range]
%    .v_error        v-wind error (m s-1)[range]
%    .w_error        w-wind error (m s-1)[range]
%    .ws_error       wind speed error (m s-1)[range]
%    .wd_error       wind direction error (degrees)[range]
%    .u_error_instr  u-wind instrumental error (m s-1)[range]
%    .v_error_instr  v-wind instrumental error (m s-1)[range]
%    .w_error_instr  w-wind instrumental error (m s-1)[range]
%    .ws_error_instr wind speed instrumental error (m s-1)[range]
%    .wd_error_instr wind direction instrumental error (m s-1)[range]
%    .R_sqred        test for horizontal homogeneity (unitless)[range]
%    .RMSE           test for horizontal homogeneity (m s-1)[range]
%    .CN             Collinearity diagnostics (unitless)[range]
%
% 2017-09-10
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

p.fit_error = true;
if ~isempty(varargin)
    p = parsePropertyValuePairs(p, varargin);
end

% Check inputs
if ~isstruct(S_in) || ~isfield(S_in,'time') || ~isfield(S_in,'range') ||...
        ~isfield(S_in,'azimuth') || ~isfield(S_in,'elevation') || ...
        ~isfield(S_in,'velocity') || ~isfield(S_in,'snr') || ~isfield(S_in,'velocity_error')
    error(sprintf(['The input must be a struct variable with followind'...
        ' fields:\n''time'', ''range'', ''azimuth'', ''elevation'',' ...
        ' ''velocity'', ''snr'', and ''velocity_error''.']))
end
if ~isnumeric(S_in.time) || isscalar(S_in.time) || ...
        any(~isfinite(S_in.time)) || any(isnan(S_in.time)) || ...
        size(S_in.time,2)>1 || ~isreal(S_in.time) || ...
        any(diff(S_in.time)<0) || any(S_in.time<0)
    error(sprintf(['The field ''time'' in ''S_in'' has to be a' ...
        ' non-negative real finite numerical vector [n 1]\nwith' ...
        ' increasing values and cannot contain any NaNs.']))
end
if ~isnumeric(S_in.range) || isscalar(S_in.range) || ...
        any(~isfinite(S_in.range)) || any(isnan(S_in.range)) || ...
        size(S_in.range,2)>1 || ~isreal(S_in.range) || ...
        any(diff(S_in.range)<0) || any(S_in.range<0)
    error(sprintf(['The field ''range'' in ''S_in'' has to be a' ...
        ' non-negative real finite numerical vector [n 1]\nwith' ...
        ' increasing values and cannot contain any NaNs.']))
end
if ~isnumeric(S_in.azimuth) || isscalar(S_in.azimuth) || ...
        any(~isfinite(S_in.azimuth(:))) || ...
        any(isnan(S_in.azimuth(:))) || size(S_in.azimuth,2)>1 || ...
        ~isreal(S_in.azimuth) || any(S_in.azimuth(:)<0) ||...
        any(S_in.azimuth(:)>360) || length(S_in.azimuth)~=length(S_in.time)
    error(sprintf(['The field ''azimuth'' in ''S_in'' has to be a real' ...
        ' finite numerical vector [time 1]\nwith all values' ...
        ' >= 0 and <= 360, and it cannot contain any NaNs.']))
end
if ~isnumeric(S_in.elevation) || isscalar(S_in.elevation) || ...
        any(~isfinite(S_in.elevation(:))) || ...
        any(isnan(S_in.elevation(:))) || size(S_in.elevation,2)>1 || ...
        ~isreal(S_in.elevation) ||  ...
        any(S_in.elevation(:)>90) || ...
        length(S_in.elevation)~=length(S_in.time)
    error(sprintf(['The field ''elevation'' in ''S_in'' has to be a' ...
        ' real finite numerical vector with all' ...
        ' values >= 0 and <= 90, and it cannot contain any NaNs.']))
end
if ~isnumeric(S_in.velocity) || isscalar(S_in.velocity) || ...
        isvector(S_in.velocity) || ~isreal(S_in.velocity) || ...
        size(S_in.velocity,1)~=length(S_in.time) || ...
        size(S_in.velocity,2)~=length(S_in.range)
    error(sprintf(['The field ''velocity'' in ''S_in'' has to be a' ...
        ' real finite numerical [time range] dimensional matrix.']))
end

if ~isnumeric(S_in.velocity_raw) || isscalar(S_in.velocity_raw) || ...
        isvector(S_in.velocity_raw) || ~isreal(S_in.velocity_raw) || ...
        size(S_in.velocity_raw,1)~=length(S_in.time) || ...
        size(S_in.velocity_raw,2)~=length(S_in.range)
    error(sprintf(['The field ''velocity_raw'' in ''S_in'' has to be a' ...
        ' real finite numerical [time range]\ndimensional matrix.']))
end

if ~isnumeric(S_in.snr) || isscalar(S_in.snr) || ...
        isvector(S_in.snr) || ~isreal(S_in.snr) || ...
        size(S_in.snr,1)~=length(S_in.time) || ...
        size(S_in.snr,2)~=length(S_in.range)
    error(sprintf(['The field ''snr'' in ''S_in'' has to be a real' ...
        ' finite numerical [time range]\ndimensional matrix.']))
end
if ~isnumeric(S_in.velocity_error) || isscalar(S_in.velocity_error) || ...
        isvector(S_in.velocity_error) || ~isreal(S_in.velocity_error) || ...
        size(S_in.velocity_error,1)~=length(S_in.time) || ...
        size(S_in.velocity_error,2)~=length(S_in.range)
    error(sprintf(['The field ''velocity_error'' in ''S_in'' has to be a' ...
        ' real finite numerical [time range]\ndimensional matrix.']))
end
% If infinite values, convert to nans
condnan = (~isfinite(S_in.snr) | ~isfinite(S_in.velocity)) | isnan(S_in.snr);
S_in.snr(condnan) = nan;
S_in.velocity(condnan) = nan;

% Initialize
u = nan(size(S_in.velocity,2),1);
v = nan(size(S_in.velocity,2),1);
w = nan(size(S_in.velocity,2),1);
u_raw = nan(size(S_in.velocity,2),1);
v_raw = nan(size(S_in.velocity,2),1);
w_raw = nan(size(S_in.velocity,2),1);
ws = nan(size(S_in.velocity,2),1);
wd = nan(size(S_in.velocity,2),1);
u_error = nan(size(S_in.velocity,2),1);
v_error = nan(size(S_in.velocity,2),1);
w_error = nan(size(S_in.velocity,2),1);
ws_error = nan(size(S_in.velocity,2),1);
wd_error = nan(size(S_in.velocity,2),1);
u_error_instr = nan(size(S_in.velocity,2),1);
v_error_instr = nan(size(S_in.velocity,2),1);
w_error_instr = nan(size(S_in.velocity,2),1);
ws_error_instr = nan(size(S_in.velocity,2),1);
wd_error_instr = nan(size(S_in.velocity,2),1);
R_squared = nan(size(S_in.velocity,2),1);
RMSE = nan(size(S_in.velocity,2),1);
CN = nan(size(S_in.velocity,2),1);

%%--- Calculate wind components ---%%
% Solve overdetermined linear system A * vr = V_r
A = [sind(S_in.azimuth(:)) .* sind(90-S_in.elevation(:)),...
    cosd(S_in.azimuth(:)) .* sind(90-S_in.elevation(:)),...
    cosd(90-S_in.elevation(:))];

[U,D,V] = svd(A); % A = U*S*transpose(V) to check

for i = 1:size(S_in.velocity,2) % loop over range gates
% Only perform fitting if enough data
    if numel(S_in.velocity(:,i)) - sum(isnan(S_in.velocity(:,i))) < 6
    % if i<4 %sum(isnan(S_in.velocity(:,i))) > numel(S_in.velocity(:,i))*.25
        u(i) = nan;
        v(i) = nan;
        w(i) = nan;
        u_raw(i) = nan;
        v_raw(i) = nan;
        w_raw(i) = nan;
        ws(i) = nan;
        wd(i) = nan;
        u_error_instr(i) = nan;
        v_error_instr(i) = nan;
        w_error_instr(i) = nan;
        ws_error_instr(i) = nan;
        wd_error_instr(i) = nan;
        R_squared(i) = nan;
        RMSE(i) = nan;
        CN(i) = nan;
        u_error(i) = nan;
        v_error(i) = nan;
        w_error(i) = nan;
        ws_error(i) = nan;
        wd_error(i) = nan;
    else
        % Wind components

        % Require finite data only - an issue for Leosphere systems with thresholding applied 
        if any(isnan(S_in.velocity(:,i)))
          % Redo SVD
          index = find(isfinite(S_in.velocity(:,i)));
          az = S_in.azimuth(index);
          el = S_in.elevation(index);
          vr = S_in.velocity(index,i);
          vr_raw = S_in.velocity_raw(index,i);
          [nU,nD,nV] = svd([sind(az(:)).*sind(90-el(:)), cosd(az(:)).*sind(90-el(:)), cosd(90-el(:))]);
          wind_components = nV * pinv(nD) * transpose(nU) * vr;
          wind_components_raw = nV * pinv(nD) * transpose(nU) * vr_raw;
        else
          wind_components = V * pinv(D) * transpose(U) * S_in.velocity(:,i);
          wind_components_raw = V * pinv(D) * transpose(U) * S_in.velocity_raw(:,i);
        end

        u(i) = wind_components(1);
        v(i) = wind_components(2);
        w(i) = wind_components(3);
        u_raw(i) = wind_components_raw(1);
        v_raw(i) = wind_components_raw(2);
        w_raw(i) = wind_components_raw(3);
        % Wind speed
        ws(i) = sqrt(u(i).^2 + v(i).^2);
        % Wind direction
        r2d = 45.0/atan(1.0); % conversion factor
        wd(i) = atan2(u(i), v(i)) * r2d + 180;
        
        %%--- Error propagation ---%        
        % variance-covariance matrices
        C_vr_vr = diag(S_in.velocity_error(:,i));
        C_v_v = pinv(A) * C_vr_vr * transpose(pinv(A));
        
        % Instrumental random noise error for u,v,w
        wind_comp_error = sqrt(diag(C_v_v));
        u_error_instr(i) = wind_comp_error(1);
        v_error_instr(i) = wind_comp_error(2);
        w_error_instr(i) = wind_comp_error(3);
        
        % Instrumental random noise error for wind speed and direction
        ws_error_instr(i) = sqrt((u(i)*u_error_instr(i))^2 + ...
            (v(i)*v_error_instr(i))^2) / sqrt(u(i)^2+v(i)^2);
        wd_error_instr(i) = sqrt((u(i)*v_error_instr(i))^2 + ...
            (v(i)*u_error_instr(i))^2) / (sqrt(u(i)^2+v(i)^2))^2;
        
        % Calculate sine wave fit and goodness-of-fit values
        if p.fit_error % true
            [~,R_squared(i),RMSE(i)] = calculateSinusoidalFit(...
                S_in.azimuth,S_in.velocity(:,i));
        else % false
            R_squared(i) = nan;
            RMSE(i) = nan;
        end
        
        % Calculate condition number
        s_1 = (transpose(A(:,1))*A(:,1))^(-1/2);
        s_2 = (transpose(A(:,2))*A(:,2))^(-1/2);
        s_3 = (transpose(A(:,3))*A(:,3))^(-1/2);
        S = diag([s_1,s_2,s_3]);
        Z = A * S;

        if any(isnan(Z(:))) | any(not(isfinite(Z(:))))
            CN(i) = nan;
        else
            CN(i) = max(svd(Z))/min(svd(Z));
        end

        
        % Radial velocity error is comprised of instrumental noise
        % and turbulent contribution since a radial velocity
        % measurement is averaged over n number of beams.
        %
        % Here the total radial velocity uncertainty is unknown,
        % and more specifically the turbulent contribution of the
        % uncertainty is unknown, and thus sigma_r is set to 1.
        % Newsom et al. (2017)
	sigma_r = 1;
        tmp_C_all = zeros(3,3);
        tmp_PSI_all = 0;
        for k = 1:length(S_in.azimuth)
            r_k = [sind((S_in.azimuth(k))).*cosd(S_in.elevation(k)),...
                cosd((S_in.azimuth(k))).*cosd(S_in.elevation(k)),...
                sind(S_in.elevation(k))];
            tmp_C = ((transpose(r_k)*r_k)/(sigma_r.^2));
            tmp_C_all = tmp_C_all + tmp_C;
            tmp_PSI = ((transpose(wind_components_raw) * transpose(r_k) - ...
                S_in.velocity_raw(k,i)).^2) / (sigma_r.^2);
            tmp_PSI_all = tmp_PSI_all + tmp_PSI;
        end
        C_newsom = inv(tmp_C_all);
        [~, mID] = lastwarn; % 'Matrix is singular to working precision' warning turned off
        if ~isempty(mID)
          warning('off',mID)
        end

        PSI = sqrt(tmp_PSI_all);
        u_error(i) = PSI * sqrt(C_newsom(1,1) / (length(S_in.azimuth) - 3));
        v_error(i) = PSI * sqrt(C_newsom(2,2) / (length(S_in.azimuth) - 3));
        w_error(i) = PSI * sqrt(C_newsom(3,3) / (length(S_in.azimuth) - 3));
        
        % % % Note: if the total uncertainty of radial velocity would be
        % % % known the wind component errors could be calculated by:
        % % % Sigma_uvw = sqrt(diag(C_newsom));
        % % % u_error(i) = Sigma_uvw(1);
        % % % v_error(i) = Sigma_uvw(2);
        % % % w_error(i) = Sigma_uvw(3);
        
        % Wind speed and direction error as given by Newsom et al. (2017)
        ws_error(i) = sqrt((u_raw(i) * u_error(i))^2 + (v_raw(i) * v_error(i))^2) / sqrt(u_raw(i)^2 + v_raw(i)^2);
        wd_error(i) = sqrt((u_raw(i) * v_error(i))^2 + (v_raw(i) * u_error(i))^2) / (sqrt(u_raw(i)^2 + v_raw(i)^2))^2;
    end
end

% Assign outputs
S_out.u = u(:);
S_out.v = v(:);
S_out.w = w(:);
S_out.ws = ws(:);
S_out.wd = wd(:);
S_out.u_error = u_error(:);
S_out.v_error = v_error(:);
S_out.w_error = w_error(:);
S_out.ws_error = ws_error(:);
S_out.wd_error = wd_error(:);
S_out.u_error_instr = u_error_instr(:);
S_out.v_error_instr = v_error_instr(:);
S_out.w_error_instr = w_error_instr(:);
S_out.ws_error_instr = ws_error_instr(:);
S_out.wd_error_instr = wd_error_instr(:);
S_out.R_squared = R_squared(:);
S_out.RMSE = RMSE(:);
S_out.CN = CN(:);







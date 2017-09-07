function [u,v,w,ws,wd,sigma_u,sigma_v,sigma_w,sigma_ws,sigma_wd,R_sqred,CN] = ...
    ppi2wind_paeshcke(site,DATE,vr,azi,ele,snr)
%PPI2WIND_PAESCHKE calculates u,v,w wind component profiles from radial velocities
% by using singular value decomposition (Paschke et al., 2015)
%
% Inputs (units) [dimensions]:
%   vr                  radial Doppler velocities (m s-1) [time,range]
%   azi                 azimuth (degrees) [time]
%   ele                 elevation (degrees) [time]
%   snr                 signal-to-noise ratio [time,range]
%
%
% Outputs (units) [dimensions]:
%   u                   u-wind profile (m s-1) [range]
%   v                   v-wind profile (m s-1) [range]
%   w                   w-wind profile (m s-1) [range]
%   ws                  wind speed profile (m s-1) [range]
%   wd                  wind direction profile (degrees) [range]
%   sigma_u             u-wind error profile (m s-1) [range]
%   sigma_v             v-wind error profile (m s-1) [range]
%   sigma_w             w-wind error profile (m s-1) [range]
%   R_squared           Goodness-of-fit for ws & wd (unitless) [range]
%   sigma_ws            wind speed error profile (m s-1) [range]
%   sigma_wd            wind direction profile (degrees) [range]
%   CN                  Condition number for ws & wd (unitless) [range]
%
% 2017-09-07
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen@helsinki.fi

C = getconfig(site,DATE);

% Initialize
u = nan(size(vr,2),1);
v = nan(size(vr,2),1);
w = nan(size(vr,2),1);
ws = nan(size(vr,2),1);
wd = nan(size(vr,2),1);
sigma_u = nan(size(vr,2),1);
sigma_v = nan(size(vr,2),1);
sigma_w = nan(size(vr,2),1);
sigma_ws = nan(size(vr,2),1);
sigma_wd = nan(size(vr,2),1);
R_sqred = nan(size(vr,2),1);
CN = nan(size(vr,2),1);

%%--- Calculate wind components ---%%
% Solve overdetermined linear system A * vr = V_r
A = [sind(azi(:)) .* sind(ele(:)),...
    cosd(azi(:)) .* sind(ele(:)),...
    cosd(ele(:))];

[U,D,V] = svd(A); % A = U*S*V' to check

% Filter background
if isfield(C,'SNR_threshold')
    vr(snr < C.snr_threshold) = nan;
else
    vr(snr < 1.01) = nan;    
end

for i = 1:size(vr,2) % loop over range gates
    if sum(isnan(vr(:,i))) > numel(vr(:,i))*.25
        u(i) = nan;
        v(i) = nan;
        w(i) = nan;
        ws(i) = nan;
        wd(i) = nan;
        sigma_u(i) = nan;
        sigma_v(i) = nan;
        sigma_w(i) = nan;
        sigma_ws(i) = nan;
        sigma_wd(i) = nan;
        R_sqred(i) = nan;
        CN(i) = nan;
    else
        wind_components = V * pinv(D) * U' * vr(:,i);
        u(i) = wind_components(1);
        v(i) = wind_components(2);
        w(i) = wind_components(3);
        
        % Wind speed 
        ws(i) = sqrt(u(i).^2 + v(i).^2);
        % Wind direction 
        r2d = 45.0/atan(1.0); % conversion factor
        wd(i) = atan2(u(i), v(i)) * r2d + 180;
        
        %%--- Error propagation ---%
        % Standard deviation of velocity estimate (Rye and Hardesty, 1997)
        obs_signal = abs(snr(:,i)-1);
        SNR_theory = 10^-5:0.01:2;
        deltav = 2;
        Np = C.num_samples_gate .* C.num_pulses_m1 .* SNR_theory';
        ff = deltav ./ (C.Nyquist .* 2);
        alpha_speckle = SNR_theory'./(sqrt(2.*pi) .* ff);
        v_err_approx = sqrt(sqrt(8).*ff.^2./(Np.*alpha_speckle).*((1+1./...
            (sqrt(2.*pi)).*alpha_speckle).^2)).*(C.Nyquist .* 2);
        sigma_Vr = interp1(SNR_theory', v_err_approx, obs_signal);
        sigma_Vr(isnan(sigma_Vr)|~isfinite(sigma_Vr)|...
            sigma_Vr > C.Nyquist) = C.Nyquist;
        
        % variance-covariance matrices
        C_vr_vr = diag(sigma_Vr);
        C_v_v = A \ C_vr_vr / A';
        
        % Calculate sigma's for u,v,w
        wind_comp_error = sqrt(diag(C_v_v));
        sigma_u(i) = wind_comp_error(1);
        sigma_v(i) = wind_comp_error(2);
        sigma_w(i) = wind_comp_error(3);
        
        % Calculate sigma's for wind speed and direction
        
        sigma_ws(i) = sqrt((u(i)*sigma_u(i))^2 + (v(i)*sigma_v(i))^2) / ...
            sqrt(u(i)^2+v(i)^2);
        sigma_wd(i) = sqrt((u(i)*sigma_v(i))^2 + (v(i)*sigma_u(i))^2) / ...
            (sqrt(u(i)^2+v(i)^2))^2;
        
        % Calculate sine wave fit's coefficient of determination R^2
        [~,R_sqred(i)] = calculateSinusoidalFit(azi,vr(:,i));
        
        % Calculate condition number
        s_1 = (A(:,1)'*A(:,1))^(-1/2);
        s_2 = (A(:,2)'*A(:,2))^(-1/2);
        s_3 = (A(:,3)'*A(:,3))^(-1/2);
        S = diag([s_1,s_2,s_3]);
        Z = A * S;
        CN(i) = max(svd(Z))/min(svd(Z));
    end
end

%%--- Nested functions ---%%
    function [sine_fit,R_squared] = calculateSinusoidalFit(azi,r_velo)
        %CALCULATESINUSOIDALFIT fits a sine wave to data.
        %https://se.mathworks.com/matlabcentral/answers/121579-curve-fitting-to-a-sinusoidal-function
        
        [~,iazi] = sort(azi); % sort w.r.t. azimuth angle (w.r.t. time orginally)
        x = azi(iazi)*pi/180; % degrees to radians
        y = r_velo(iazi)-nanmean(r_velo); % shift by mean
        
        yu = max(y);
        yl = min(y);
        % Range of ?y?
        yr = (yu-yl);
        yz = y-yu+(yr/2);
        % Find zero-crossings
        zx = x(yz .* circshift(yz,[-1 1]) <= 0);
        % Estimate period
        per = 2*mean(diff(zx));
        % Estimate offset
        ym = mean(y);
        
        % Function to fit
        fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);
        % Least-Squares cost function
        fcn = @(b) sum((fit(b,x) - y).^2);
        % Minimise Least-Squares
        s = my_fminsearch(fcn, [yr;  per;  -1;  ym]);
        
        xp = linspace(min(x),max(x),length(x));
        
        % figure; hold on;
        % plot(x*(180/pi),y+nanmean(r_velo(ia)),'b')
        % plot(xp*(180/pi),fit(s,xp)+nanmean(r_velo(ia)), 'r')
        % grid
        
        sine_fit = fit(s,xp)+nanmean(r_velo); % calc fit and shift back
        [~,itime] = sort(iazi);
        sine_fit = sine_fit(itime); % order back w.r.t time
        R_squared = 1 - sum((r_velo(iazi) - ...
            (fit(s,xp)+nanmean(r_velo(iazi)))'  ).^2) / ...
            sum((r_velo(iazi) - sum(r_velo(iazi))).^2);
    end

    function [x,fval] = my_fminsearch(funfcn,x)
        
        maxfun = 800;
        maxiter = 800;
        tolf = 1.0000e-4;
        tolx = 1.0000e-4;
        varargin = {};
        n = numel(x);
        
        % Initialize parameters
        rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
        onesn = ones(1,n);
        two2np1 = 2:n+1;
        one2n = 1:n;
        
        % Set up a simplex near the initial guess.
        xin = x(:); % Force xin to be a column vector
        vee = zeros(n,n+1); fv = zeros(1,n+1);
        vee(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
        x(:) = xin;    % Change x to the form expected by funfcn
        fv(:,1) = funfcn(x,varargin{:});
        func_evals = 1;
        itercount = 0;
        how = '';
        % Initial simplex setup continues later
        
        % Continue setting up the initial simplex.
        % Following improvement suggested by L.Pfeffer at Stanford
        usual_delta = 0.05;             % 5 percent deltas for non-zero terms
        zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
        for j = 1:n
            y = xin;
            if y(j) ~= 0
                y(j) = (1 + usual_delta)*y(j);
            else
                y(j) = zero_term_delta;
            end
            vee(:,j+1) = y;
            x(:) = y; f = funfcn(x,varargin{:});
            fv(1,j+1) = f;
        end
        
        % sort so vee(1,:) has the lowest function value
        [fv,j] = sort(fv);
        vee = vee(:,j);
        
        how = 'initial simplex';
        itercount = itercount + 1;
        func_evals = n+1;
        
        
        % Main algorithm: iterate until
        % (a) the maximum coordinate difference between the current best point and the
        % other points in the simplex is less than or equal to TolX. Specifically,
        % until max(||v2-v1||,||v3-v1||,...,||vee(n+1)-v1||) <= TolX,
        % where ||.|| is the infinity-norm, and v1 holds the
        % vertex with the current lowest value; AND
        % (b) the corresponding difference in function values is less than or equal
        % to TolFun. (Cannot use OR instead of AND.)
        % The iteration stops if the maximum number of iterations or function evaluations
        % are exceeded
        while func_evals < maxfun && itercount < maxiter
            if max(abs(fv(1)-fv(two2np1))) <= max(tolf,10*eps(fv(1))) && ...
                    max(max(abs(vee(:,two2np1)-vee(:,onesn)))) <= max(tolx,10*eps(max(vee(:,1))))
                break
            end
            
            % Compute the reflection point
            
            % xbar = average of the n (NOT n+1) best points
            xbar = sum(vee(:,one2n), 2)/n;
            xr = (1 + rho)*xbar - rho*vee(:,end);
            x(:) = xr; fxr = funfcn(x,varargin{:});
            func_evals = func_evals+1;
            if fxr < fv(:,1)
                % Calculate the expansion point
                xe = (1 + rho*chi)*xbar - rho*chi*vee(:,end);
                x(:) = xe; fxe = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                if fxe < fxr
                    vee(:,end) = xe;
                    fv(:,end) = fxe;
                    how = 'expand';
                else
                    vee(:,end) = xr;
                    fv(:,end) = fxr;
                    how = 'reflect';
                end
            else % fv(:,1) <= fxr
                if fxr < fv(:,n)
                    vee(:,end) = xr;
                    fv(:,end) = fxr;
                    how = 'reflect';
                else % fxr >= fv(:,n)
                    % Perform contraction
                    if fxr < fv(:,end)
                        % Perform an outside contraction
                        xc = (1 + psi*rho)*xbar - psi*rho*vee(:,end);
                        x(:) = xc; fxc = funfcn(x,varargin{:});
                        func_evals = func_evals+1;
                        
                        if fxc <= fxr
                            vee(:,end) = xc;
                            fv(:,end) = fxc;
                            how = 'contract outside';
                        else
                            % perform a shrink
                            how = 'shrink';
                        end
                    else
                        % Perform an inside contraction
                        xcc = (1-psi)*xbar + psi*vee(:,end);
                        x(:) = xcc; fxcc = funfcn(x,varargin{:});
                        func_evals = func_evals+1;
                        
                        if fxcc < fv(:,end)
                            vee(:,end) = xcc;
                            fv(:,end) = fxcc;
                            how = 'contract inside';
                        else
                            % perform a shrink
                            how = 'shrink';
                        end
                    end
                    if strcmp(how,'shrink')
                        for j=two2np1
                            vee(:,j)=vee(:,1)+sigma*(vee(:,j) - vee(:,1));
                            x(:) = vee(:,j); fv(:,j) = funfcn(x,varargin{:});
                        end
                        func_evals = func_evals + n;
                    end
                end
            end
            [fv,j] = sort(fv);
            vee = vee(:,j);
            itercount = itercount + 1;
        end   % while
        
        x(:) = vee(:,1);
        fval = fv(:,1);
    end


end


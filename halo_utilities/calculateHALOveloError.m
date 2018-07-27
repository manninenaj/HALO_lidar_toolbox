function [v_error, v_error_att] = calculateHALOveloError(site,DATE,measmode,typeof,signal)
%calculateHALOveloError calculates the standard deviation of an individual 
% velocity estimate.
% 
% Inputs:
% - site            string, name of the site, e.g. 'kuopio'
% - DATE            scalar, numerical date, e.g. 20171231
% - snr             matrix, signal-to-noise ratio, no. of rows must match 
%                   the length of time
%
% Ouputs:
% - v_error         Standard deviation of the velocity estimate, calculated
%                   from the SNR variable following Rye and Hardesty (1997)
% - v_error_att     v_error attirubtes

% Check inputs
if nargin < 5
    error(['''site'', ''DATE'', ''measmode'', ''typeof'', and ''signal'''...
        ' are required inputs!'])
end
if ~ischar(site)
    error('The 1st input ''site'' must be a string.')
end
if ~isnumeric(DATE) || length(num2str(DATE))~=8
    error(['The 2nd input ''DATE'' must be a numerical date in' ...
        ' YYYYMMDD format.'])
end
if ~ischar(measmode) || ~any(strcmp(measmode,{'stare','vad','dbs','rhi','custom'}))
    error(sprintf(['The 3rd input ''measmode'' must be a string and can be:\n'...
        '''stare'',''vad'',''rhi'',''dbs'',''custom''.']))
end
if isscalar(signal) || ~isnumeric(signal)
    error('The 5th input ''signal'' must a numerical matrix.')
end

% Get parameters
C = getconfig(site,DATE);
abc = [measmode '_' typeof];

theory.PRF = C.prf; % pulse repetition frequency
theory.Nyquist = C.nyquist; % Nyquist velocity
theory.B = theory.Nyquist .* 2;
theory.M = C.(['num_samples_gate_' abc]); % points per range gate
theory.npulses = C.(['num_pulses_m1_' abc]);
obs_signal = abs(signal-1);
% theory.SNR = transpose([10.^[-5:0.01:2]]);
theory.SNR = transpose([10.^[-5:0.01:2.5]]);
theory.deltav = [1 1.5 2]; % typical signal spectral width - best value
                           % 1.5, or 2?
theory.ff = theory.deltav./theory.B;
theory.Np = theory.M .* theory.npulses .* theory.SNR;
for ii = 1:length(theory.deltav)
  % alpha = ratio of photon count to speckle count
  % alpha_speckle = mySNR.*6;
  theory.alpha_speckle(:,ii) = theory.SNR./(sqrt(2.*pi).*theory.ff(ii));
  theory.v_err_approx(:,ii) = sqrt(sqrt(8).*theory.ff(ii).^2./...
      (theory.Np .*theory.alpha_speckle(:,ii)).*((1+1./(sqrt(2.*pi)).*...
      theory.alpha_speckle(:,ii)).^2)).*theory.B;
  theory.direct_detection(:,ii) = theory.deltav(ii)./sqrt(theory.Np);
  theory.v_err_pearson(:,ii) = 2.*sqrt(sqrt(pi)./...
      theory.alpha_speckle(:,ii)).*(theory.deltav(ii)./...
      sqrt(theory.Np)).*(1+(1./(2.*pi)).*theory.alpha_speckle(:,ii));
  % full treatment of gaussian pulse requires numerical integration
  % The limits should be -0.5 to 0.5. Integrate from 0 to 0.5 and
  % multiply by 2. This provides the Cramer-Rao lower bound (CRLB)
  x = 0:0.001:0.5;
  a = 10.^[-4:0.1:4];  
  g1 = zeros(size(a));
  f = zeros(size(a));
  for ialpha = 1:length(a)
    g1(ialpha) = a(ialpha)./sqrt(2.*pi).*2.*sum((x.^2.*exp(-(x.^2)))./...
        (1+a(ialpha).*(exp(-(x.^2)./2).^2))).*0.001;
    f(ialpha) = 2.*sum(((x./theory.ff(ii)).^2)./((1 + 1./(a(ialpha).*...
        (exp(-(x.^2)./(2.*theory.ff(ii).^2))))).^2)).*0.001;
  end
  theory.f(:,ii) = interp1(a,f,theory.alpha_speckle(:,ii));
  theory.v_err_crlb(:,ii) = sqrt(theory.ff(ii).^2./theory.npulses./...
      theory.M./theory.f(:,ii)).*theory.B;
end
v_error = interp1(theory.SNR, theory.v_err_crlb(:,3), obs_signal);
v_error(isnan(v_error)|~isfinite(v_error)|v_error > theory.Nyquist) = ...
    theory.Nyquist;
v_error_att = create_attributes({'time','range'},'Std of velocity',...
    {'ms-1','m s<sup>-1</sup>'}, C.missing_value, ['This variable is' ...
    ' the standard deviation of the velocity estimate, calculated' ...
    ' from the SNR variable following Rye and Hardesty (1997). Note' ...
    ' that this value is capped at the Nyquist velocity.']);

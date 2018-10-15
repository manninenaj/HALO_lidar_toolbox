function [pr2,pr2_0] = correct_focus(focus, data)

c_light =  299792458; 
hnu = 1.28e-19;
alpha = 0.01;

if ~isfield(data, 'lens_diameter')
  data.lens_diameter = 0.06;
end

if ~isfield(data, 'wavelength')
  data.wavelength = 1.5e-6;
end

if ~isfield(data, 'beam_energy')
  data.beam_energy = 0.00001; % 1e-5
end

het_area = pi .* ((0.7 * data.lens_diameter ./ 2).^2);

if focus < 0
  nt = het_area ./ data.wavelength ./ data.range;
else
  nt = het_area ./ data.wavelength .* (1 ./ focus  - 1 ./ data.range);
end

effective_area = het_area ./ (1 + nt.^2);
T = 10.^(-1 .* alpha .* data.range / 5000);
het_cnr = 0.7 .* 0.4 .* 0.6 .* effective_area .* c_light .* data.beam_energy .* T ./ (2 .* data.range.^2 .* hnu .* data.bandwidth);

pr2 = (data.signal - 1) ./ (ones(length(data.time),1) * het_cnr' );
pr2_0 = (data.signal0 - 1) ./ (ones(length(data.time),1) * het_cnr' );

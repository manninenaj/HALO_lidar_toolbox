function [sine_fit,R_squared,RMSE] = calculateSinusoidalFit(azi,vr)
%calculateSinusoidalFit fits a sine wave to data.

if all(isnan(vr))

  sine_fit = nan(size(azi));
  R_squared = nan;
  RMSE = nan(size(azi));
  RMSE = RMSE(:); R_squared = R_squared(:); sine_fit = sine_fit(:); 

else

  [~,iazi] = sort(azi); % sort w.r.t. azi (w.r.t. time orginally)
  x = azi(iazi)*pi/180; % degrees to radians
  y = vr(iazi)-nanmean(vr); % shift by mean

  yu = nanmax(y);
  yl = nanmin(y);

  % Range of ?y?
  yr = (yu-yl);
  %yz = y-yu+(yr/2);
  %%% Find zero-crossings
  %%% zx = x(yz .* circshift(yz,[-1 1]) <= 0);
  % Estimate period
  %%% per = 2*nanmean(diff(zx));
  per = 2.*pi;
  % Estimate offset
  ym = nanmean(y);

  % Function to fit
  func_to_fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);
  % Least-Squares cost function
  fcn = @(b) nansum((func_to_fit(b,x) - y).^2);
  % Minimise Least-Squares
  s = funcMinimizer(fcn, [yr;  per;  -1;  ym]);
  xp = linspace(nanmin(x),nanmax(x),length(x));

  % calc fit and shift back
  sine_fit = func_to_fit(s,xp) + nanmean(vr);
  [~,itime] = sort(iazi);
  sine_fit = sine_fit(itime); % order back w.r.t time

  % Goodness-of-fit parameters
  R_squared = 1 - nansum((vr(iazi) - (transpose(func_to_fit(s,xp)) + nanmean(vr(iazi)))  ).^2) ./ nansum((vr(iazi) - nansum(vr(iazi))).^2);

  RMSE = sqrt(nansum(sine_fit(:) - vr(:)).^2 / length(vr));
end



function [xp,yp,zp] = pcolor2polar(time,range,azimuth,data)

  % Original data (in cartesian space)
  x = time;
  y = range;
  z = data;
  % Convert to polar
  r = range;
  th = azimuth;
  for i = 1:length(r)
    for j = 1:length(th)
      xp(i,j) = r(i).*cosd(th(j));
      yp(i,j) = r(i).*sind(th(j));
    end
  end
 whos 
  zp = interp2(x,y,z,xp',yp');
  %[thpt, rpt] = cart2pol(xpt, ypt);
  
  %% Plot in polar
  %ax(2) = subplot(1,2,2);
  %pcolor(xp,yp,zp);
  %shading flat;
  %axis equal;
  %ax(3) = axes('position', get(ax(2), 'position'));
  %hp = mmpolar(thpt, rpt, 'ko', 'RLimit', [0 50], 'backgroundcolor', 'none');
  %set(ax(2), 'visible', 'off');

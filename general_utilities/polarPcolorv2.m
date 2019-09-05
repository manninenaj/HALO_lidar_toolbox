function polarPcolorv2(time,range,azimuth,data)

  % Original data (in cartesian space)
  x = time;
  y = range;
  z = data;

  tmp1 = rand(3,1)*2*pi;
  tmp2 = rand(3,1)*50;
  % Convert to polar
  r = range;
  th = azimuth;
  xp = r.*cosd(th);
  yp = r.*sind(th);
  zp = interp2(x,y,z,xp,yp);
  [thpt, rpt] = cart2pol(xpt, ypt);
  
  % Plot in polar
  ax(2) = subplot(1,2,2);
  pcolor(xp,yp,zp);
  shading flat;
  axis equal;
  ax(3) = axes('position', get(ax(2), 'position'));
  hp = mmpolar(thpt, rpt, 'ko', 'RLimit', [0 50], 'backgroundcolor', 'none');
  set(ax(2), 'visible', 'off');

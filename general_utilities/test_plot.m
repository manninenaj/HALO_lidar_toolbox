% Original data (in cartesian space)
  [x,y] = meshgrid(linspace(-50,50,200));
z = peaks(200);
tmp1 = rand(3,1)*2*pi;
tmp2 = rand(3,1)*50;
xpt = tmp2 .* cos(tmp1);
ypt = tmp2 .* sin(tmp1);
% Plot in cartesian
ax(1) = subplot(1,2,1);
pcolor(x,y,z);
shading flat;
hold on;
plot(xpt, ypt, 'ko');
axis equal tight;
% Convert to polar
[r,th] = meshgrid(linspace(0,50,100), linspace(0, 2*pi, 360));
xp = r.*cos(th);
yp = r.*sin(th);
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

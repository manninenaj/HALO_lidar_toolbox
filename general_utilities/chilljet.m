function [color_map] = chilljet(varargin)

if length(varargin) > 1
    error 'Too many inputs!'
end

% Create colormap base
color_map_base =  [...
128  0  128;...
256 0   0;...
256 165 0;...
256 256 0;...
0   256 0;...
32  178 170;...
0   0   256;...
100 149 237;...
230 230 250];

% Expand to higher reso color map
color_map = nan((size(color_map_base,1)-1)*90, size(color_map_base,2));
for iR = 1:size(color_map_base,1)-1 % rows index
    for iC = 1:size(color_map_base,2) % col index
        if iR ~= 1 % Adjust indeces
            iR_s = (iR-1) * 100;
            iR_e = (iR-1) * 100 + 100;
            n = 101;
        else  % Adjust indeces
            iR_s = iR;
            iR_e = iR + 99;
            n = 100;
        end
    
        color_map(iR_s:iR_e,iC) = linspace(color_map_base(iR,  iC),...
            color_map_base(iR+1,iC),n);

    end
end

% Convert from 0 to 1
color_map = flipud(color_map ./ 256);

if length(varargin) == 1
    tmp = nan(varargin{1},3);
    tmp(:,1) = interp1(1:length(color_map(:,1)),color_map(:,1),linspace(1,length(color_map(:,1)),varargin{1}));
    tmp(:,2) = interp1(1:length(color_map(:,2)),color_map(:,2),linspace(1,length(color_map(:,2)),varargin{1}));
    tmp(:,3) = interp1(1:length(color_map(:,3)),color_map(:,3),linspace(1,length(color_map(:,3)),varargin{1}));
    color_map = tmp;
end
end


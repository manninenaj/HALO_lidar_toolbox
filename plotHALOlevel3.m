function plotHALOlevel3(site,DATES,processing_level,observation_type,varargin)

p.sub_type = nan;
p.ylabel = '';
p.xlabel = 'Time UTC (hrs)';
p.masking = 1; % SNR+1
p.ylim = [0 nan]; % km
p.ystep = 2; % km
p.azilim = [0 360]; % degrees
p.azistep = 60; % degrees
p.cmap = cmap_darkviolet_to_brickred;
p.cmapdiv = cmocean('balance');
p.cmapwdir = colorcet('C8');
if ~isempty(varargin)
    p = parsePropertyValuePairs(p, varargin);
    p.ystep = p.ystep/1000; % km
    if p.ylim(1) ~= 0
        p.ylim(1) = p.ylim(1)/1000; % km
    end
    p.ylim(2) = p.ylim(2)/1000; % km
    if isfield(p,'sub_type')
      sub_type = p.sub_type;
    end
end

% Check inputs
list_of_processing_levels = {'original','corrected','calibrated','background','product','level3'};
list_of_observation_types = {'stare','vad','dbs','rhi','custom','co','txt','nc'}; % TODO: make txt and nc optional inputs
list_of_products = {'windvad','winddbs','epsilon','wstats','wstats4precipfilter','sigma2vad','windshear',...
    'LLJ','ABLclassification','cloud','betavelocovariance','ABLclassificationClimatology'};
if nargin < 4
    error("At least inputs 'site', 'DATE', 'processing_level', and 'observation_type'")
end
if (nargin == 4 || nargin == 5) && (strcmp(processing_level,'product') && any(strcmp(observation_type,list_of_products)) || ...
       strcmp(processing_level,'background'))
    if ~ischar(site)
        error("The 1st input (site name) must be a string.")
    end
    if ~isnumeric(DATE) || length(num2str(DATE))~=8
        error("The 2nd input (date) must be a numeric value in YYYYMMDD format.")
    end
    if ~ischar(processing_level) || ~any(strcmp(processing_level, list_of_processing_levels))
        error("The 3rd input (processing level) must be a string and one of these:\n%s", ...
           sprintf('%s,', list_of_processing_levels{:}))
    end
    if ~ischar(observation_type) || ~any(strcmp(observation_type,[list_of_observation_types, list_of_products]))
        error("The 4th input (observation type) must be a string and one of these:\n%s", ...
            sprintf("'%s','%s',", list_of_observation_types{:}, list_of_products{:}))
    end
end
if nargin < 5 && (~strcmp(processing_level,'product') && ~any(strcmp(observation_type, list_of_products)) && ...
       ~strcmp(processing_level,'background'))
    error("Input 'sub_type' is required with %s'", observation_type)
end

hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 .5 25 10];
hf.Color = 'white'; hf.Visible = 'off';
% hf = figure; hf.Units = 'Normalized'; hf.Position = [.2 .1 .4 .65]; hf.Color = 'w';
% % summer
% sp3 = subplot(223);
% bar(N_blclass_summer','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'c)'); title('JJA')
% set(gca,'XTick',2:2:24); xlabel('Time of day UTC'); ylabel('Probability')
% sp3.Position = [.075 .065 .4 .32];
% 
% % autumn
% sp4 = subplot(224);
% bar(N_blclass_autumn','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'd)'); title('SON')
% set(gca,'XTick',2:2:24); xlabel('Time of day UTC');
% sp4.Position = [.55 .065 .4 .32];
% 
% % winter
% sp1 = subplot(221);
% bar(N_blclass_winter','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'a)'); title('DJF')
% set(gca,'XTick',2:2:24); ylabel('Probability')
% sp1.Position = [.075 .48 .4 .32];
% 
% % spring
% sp2 = subplot(222);
% bar(N_blclass_spring','stacked','EdgeColor','none'); shading flat;
% colormap([blclass_red(:) blclass_green(:) blclass_blue(:)])
% axis([.5 24.5 0 1]); text(.5, 1.05, 'b)'); title('MAM')
% set(gca,'XTick',2:2:24);
% sp2.Position = [.55 .48 .4 .32];
% 
% 
% M = {'Missing data','Non-turbulent','Convective mixing','Wind shear',...
%     'Decaying / intermittent','In cloud','Cloud driven'};
% hl = legend(M); hlpos = hl.Position; hlpos(1:2) = [.4 .82];
% hl.Position = hlpos; hl.Box = 'off';
end

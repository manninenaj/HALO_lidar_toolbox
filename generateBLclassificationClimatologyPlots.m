function generateBLclassificationClimatologyPlots(file_name,varargin)

% p.height_bin_size = 100; % m
month_names = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
p.seasons = {[12,1,2],[3,4,5],[6,7,8],[9,10,11]};
p.ylabel = 'Height (km)';
p.xlabel = 'TOD UTC (hrs)';
p.ylims = [0 2.5]; % km
p.y_tick_step = 1; % km
p.cmap = cmap_darkviolet_to_brickred;
if ~isempty(varargin)
    p = parsePropertyValuePairs(p, varargin);
    p.y_tick_step = p.ystep;
    if p.ylims(1) ~= 0
        p.ylims(1) = p.ylims(1);
    end
    p.ylims(2) = p.ylims(2);
    if length([p.seasons]) ~= 12
        error("'season' parameters incorrect!")
    end
end

% Check inputs
if ~ischar(file_name)
    error('The first input must be a string.')
end

% Use datenum to accommodate leap years etc.
[data, att] = load_nc_struct(file_name);
x = data.hour;
y = data.height_agl/1000;

for i = 1:length(data.bl_classification_types)
   istart = strfind(att.bl_classification_types.definition,':')-1;
   istop = [istart(2:end) length(att.bl_classification_types.definition)];
   bl_class_names{i} = att.bl_classification_types.definition(istart(i):istop(i)-2);
end
fnames = fieldnames(data);
fnames_time = fnames(strmatch('counts_bl_classification', fnames));
for ir = 1:length(fnames_time)
    strtmp = strfind(fnames_time{ir},'_');
    tres = fnames_time{ir}(strtmp(end)+1:end);
    
    n_classes = length(data.bl_classification_types)-1;
    n_seasons = length(p.seasons);
    
    hf = figure;
    hf.Units = 'centimeters';
    hf.Position = [.5 .5 20 20];
    hf.Color = 'white';
    hf.Visible = 'off';
    
    ip = 1;
    for iclass = 1:n_classes
        for iseason = 1:n_seasons
            counts = 0;
            elements = 0;
            for imon = p.seasons{iseason}
                counts = counts + squeeze(data.(['counts_bl_classification_' tres])(iclass,imon,:,:));
                elements = elements + squeeze(data.(['number_of_elements_' tres])(imon,:,:));
            end
            
            % Calculate probability
            z = counts ./ elements;
            
            subplot(n_classes, n_seasons, ip);
            pcolor(x, y, z');
            shading flat;
            caxis([0 1]);
            colormap(p.cmap)
            axis([0 24 0 p.ylims(2)])
            set(gca,'XTick',0:6:24,'YTick',0:p.y_tick_step:p.ylims(2))
            if ip == n_seasons*iclass
                pause(.2)
                apos = get(gca,'Position');
                cb = colorbar;
                cb.Units = 'centimeters';
                cb.Position(3) = .1;
                cb.Position(1) = 18.5;
                cb.Ticks = 0:.2:1;
                cb.Label.String = 'probability';
                set(gca,'Position',apos);
            end
            if iclass == 1
                tmp = month_names(p.seasons{iseason});
                title([tmp{:}])
            end
            if iseason == 1
                ylabel(p.ylabel)
            end
            if iclass == n_classes
                xlabel(p.xlabel)
            end
            text(.75,p.ylims(2)*.85,bl_class_names{iclass},'Color','w','FontSize',6)
            ip = ip + 1;
        end
    end
end

file_name_out = strrep(file_name,'.nc','.png');
fprintf('Writing %s\n',file_name_out)
export_fig('-png','-m2','-nocrop',file_name_out)

close(hf)


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

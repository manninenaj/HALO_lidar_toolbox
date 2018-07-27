function calculateHALOgroupedDataStats(site,DATES)

% Check inputs
if nargin < 2
    error('''site'' and ''DATES'' are required inputs!')
end
if ~ischar(site)
    error('The first input ''site'' must be a string.')
end
if length(DATES)>2
    error('''DATES'' can have max. length of 2.')
elseif length(DATES)==1
    if length(num2str(DATES))~=8
        error(['The value in the second input ''DATES'' must be' ...
            ' numerical date in YYYYMMDD format.'])
    else
        DATEstart = DATES; DATEend = DATES;
    end
elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
        length(num2str(DATES(2)))~=8)
    error(['The value(s) in the second input ''DATES'' must be' ...
        ' numerical date(s) in YYYYMMDD format.'])
else
    DATEstart = DATES(1); DATEend = DATES(2);
end

% Define edges and bin centers
edges_skewn = linspace(-6,6,200);
binc_skewn = edges_skewn - median(diff(edges_skewn))/2;
binc_skewn(1) = [];
edges_tke = linspace(-9,2,200);
binc_tke = edges_tke - median(diff(edges_tke))/2;
binc_tke(1) = [];

% Initialize
N_100m_skewn = zeros(length(edges_skewn)-1,24,12); % bins-hrs-months
N_300m_skewn = zeros(length(edges_skewn)-1,24,12);
N_800m_skewn = zeros(length(edges_skewn)-1,24,12);
N_100m_tke = zeros(length(edges_tke)-1,24,12); % bins-hrs-months
N_300m_tke = zeros(length(edges_tke)-1,24,12);
N_800m_tke = zeros(length(edges_tke)-1,24,12);

tke_med_100m = nan(24,12);
tke_med_300m = nan(24,12);
tke_med_800m = nan(24,12);
skewn_med_100m = nan(24,12);
skewn_med_300m = nan(24,12);
skewn_med_800m = nan(24,12);

tke_p25_100m = nan(24,12);
tke_p25_300m = nan(24,12);
tke_p25_800m = nan(24,12);
skewn_p25_100m = nan(24,12);
skewn_p25_300m = nan(24,12);
skewn_p25_800m = nan(24,12);

tke_p75_100m = nan(24,12);
tke_p75_300m = nan(24,12);
tke_p75_800m = nan(24,12);
skewn_p75_100m = nan(24,12);
skewn_p75_300m = nan(24,12);
skewn_p75_800m = nan(24,12);

% Use datenum to accommodate leap years
for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):datenum(num2str(DATEend),'yyyymmdd')
    
    % Convert back to numerical date
    DATE = str2double(datestr(DATEi,'yyyymmdd'));
    thedate = datestr(DATEi,'yyyymmdd');
    [~,i_month] = datevec(DATEi);
  
    [dir_wstats,files_wstats] = getHALOfileList(site,DATE,'product','wstats');
    [dir_tke,files_tke] = getHALOfileList(site,DATE,'product','TKE');
    
    if isempty(files_wstats) || isempty(files_tke), continue; end
    
    wstats = load_nc_struct_silent(fullfile([dir_wstats files_wstats{1}]),...
        {'time_30min','height','radial_velocity_weighted_skewness_30min','radial_velocity_weighted_skewness_error_30min'});
    tke = load_nc_struct_silent(fullfile([dir_tke files_tke{1}]),...
	{'time_3min','height','epsilon_w_3min','epsilon_w_error_3min'});

    epsilon = tke.epsilon_w_3min;
    epsilon(tke.epsilon_w_error_3min > 1) = nan;
    skewn = wstats.radial_velocity_weighted_skewness_30min;
    skewn_e = wstats.radial_velocity_weighted_skewness_error_30min;
    skewn(real(log10(skewn_e))>-.1) = nan;

    for i_hr = 0:23
        % Extract
        skewn_100m = skewn(floor(wstats.time_30min) == i_hr,4:6);
        skewn_300m = skewn(floor(wstats.time_30min) == i_hr,10:12);
        skewn_800m = skewn(floor(wstats.time_30min) == i_hr,27:29);
        % Counts
        N_100m_skewn_tmp = hist(skewn_100m(:),binc_skewn);
        N_300m_skewn_tmp = hist(skewn_300m(:),binc_skewn);
        N_800m_skewn_tmp = hist(skewn_800m(:),binc_skewn);
        % Add up
        N_100m_skewn(:,i_hr+1,i_month) = N_100m_skewn(:,i_hr+1,i_month) + N_100m_skewn_tmp(:);
        N_300m_skewn(:,i_hr+1,i_month) = N_300m_skewn(:,i_hr+1,i_month) + N_300m_skewn_tmp(:);
        N_800m_skewn(:,i_hr+1,i_month) = N_800m_skewn(:,i_hr+1,i_month) + N_800m_skewn_tmp(:);
        % Extract
        tke_100m = real(log10(epsilon(floor(tke.time_3min) == i_hr,4:6)));
        tke_300m = real(log10(epsilon(floor(tke.time_3min) == i_hr,10:12)));
        tke_800m = real(log10(epsilon(floor(tke.time_3min) == i_hr,27:29)));
        % Counts
        N_100m_tke_tmp = hist(tke_100m(:),binc_tke);
        N_300m_tke_tmp = hist(tke_300m(:),binc_tke);
        N_800m_tke_tmp = hist(tke_800m(:),binc_tke);
        % Add up
        N_100m_tke(:,i_hr+1,i_month) = N_100m_tke(:,i_hr+1,i_month) + N_100m_tke_tmp(:);
        N_300m_tke(:,i_hr+1,i_month) = N_300m_tke(:,i_hr+1,i_month) + N_300m_tke_tmp(:);
        N_800m_tke(:,i_hr+1,i_month) = N_800m_tke(:,i_hr+1,i_month) + N_800m_tke_tmp(:);
    end
end

for i_hr2 = 1:24
    for i_m = 1:12
        % Calculate median, 25th, 75th percentile from the grouped data
        [tke_med_100m(i_hr2,i_m), tke_p25_100m(i_hr2,i_m), tke_p75_100m(i_hr2,i_m)] = ...
            groupedDataStats(N_100m_tke(:,i_hr2,i_m),binc_tke(:));
        [tke_med_300m(i_hr2,i_m), tke_p25_300m(i_hr2,i_m), tke_p75_300m(i_hr2,i_m)] = ...
            groupedDataStats(N_300m_tke(:,i_hr2,i_m),binc_tke(:));
        [tke_med_800m(i_hr2,i_m), tke_p25_800m(i_hr2,i_m), tke_p75_800m(i_hr2,i_m)] = ...
            groupedDataStats(N_800m_tke(:,i_hr2,i_m),binc_tke(:));

        [skewn_med_100m(i_hr2,i_m), skewn_p25_100m(i_hr2,i_m), skewn_p75_100m(i_hr2,i_m)] = ...
            groupedDataStats(N_100m_skewn(:,i_hr2,i_m),binc_skewn(:));
        [skewn_med_300m(i_hr2,i_m), skewn_p25_300m(i_hr2,i_m), skewn_p75_300m(i_hr2,i_m)] = ...
            groupedDataStats(N_300m_skewn(:,i_hr2,i_m),binc_skewn(:));
        [skewn_med_800m(i_hr2,i_m), skewn_p25_800m(i_hr2,i_m), skewn_p75_800m(i_hr2,i_m)] = ...
            groupedDataStats(N_800m_skewn(:,i_hr2,i_m),binc_skewn(:));
    end
end

%%

monthnames = cellstr(datestr(datenum(2017,1,1):31:datenum(2017,12,31),'mmm'));
hf = figure; hf.Units = 'Normalized'; hf.Position = [0 0 1 1]; hf.Color = 'w';
for m = 1:12
    if strcmp(site,'arm-graciosa') && m < 6
        subplot(3,4,m);
        hold on
        
        %     p1 = patch([0:23 23:-1:0], [smooth(p25_latent.(['m' num2str(m)]),3,'rlowess'); flipud(smooth(p75_latent.(['m' num2str(m)]),3,'rlowess'))], rgb('DarkGray'),'EdgeColor','none');
        %     alpha(0.4)
        ax1 = gca; % current axes
        ax1.YColor = rgb('Crimson');
        ax1.XTick = 0:3:24;
        ax1.XLim = [0 24];
        ax1.YLim = [-5 -2];
        ax1_pos = ax1.Position; % position of first axes
        ax1.Box = 'on';
        ax1.YTickLabels = [repmat('10^{',length(ax1.YTick(:)),1) num2str(ax1.YTick(:)) repmat('}',length(ax1.YTick(:)),1)];
        
        %     e1 = line((0:23),med_latent.(['m' num2str(m)]),'Parent',ax1);
        switch m
            case {1,2,3,4,5}
                e1 = line((0.5:23.5),tke_med_100m(:,m),'Parent',ax1);
            case {6,7,8,9,10}
                e1 = line((0.5:23.5),tke_med_100m(:,m),'Parent',ax1);
            case {11,12,13,14,15}
                e1 = line((0.5:23.5),tke_med_100m(:,m),'Parent',ax1);
        end
        
        set(e1,'LineStyle',':','Marker','.','MarkerSize',10,'Color',rgb('Crimson'));
        
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','right',...
            'Color','none');
        ax2.YLim = [-1 1];
        ax2.XTick = 0:3:24;
        ax2.XLim = [0 24];
        ax2.YColor = rgb('Black');
        %     p2 = patch([0:23 23:-1:0], [smooth(p25_sensible.(['m' num2str(m)]),3,'rlowess'); flipud(smooth(p75_sensible.(['m' num2str(m)]),3,'rlowess'))], rgb('DarkSalmon'),'EdgeColor','none');
        %     alpha(0.4)
        
        switch m
            case {1,2,3,4,5}
                e2 = line((0.5:23.5),skewn_med_100m(:,m),'Parent',ax2);
            case {6,7,8,9,10}
                e2 = line((0.5:23.5),skewn_med_100m(:,m),'Parent',ax2);
            case {11,12,13,14,15}
                e2 = line((0.5:23.5),skewn_med_100m(:,m),'Parent',ax2);
        end
        set(e2,'LineStyle',':','Marker','.','MarkerSize',10,'Color',rgb('Black'));
        
        grid on;
        switch m
            case {9,10,11,12}
                xlabel('Time of day (UTC)');
        end
        switch m
            case {1,5,9}
                ylabel(ax1,'TKE (m^{2} s^{-3})');
        end
        switch m
            case {4,8,12}
                ylabel(ax2,'Skewness (unitless)');
        end
        if strcmp(site,'arm-graciosa') && m == 5
            ylabel(ax2,'Skewness (unitless)');
        end
        hold off
        title(monthnames{m})
    end
end
[hl] = legend([e1 e2],{' TKE,','Skewness'},'Orientation','Horizontal','Box','On');
hl.Position = [.275 .96 .5 .02];

% export_fig -nocrop -png 20160907-20170922_arm-ascension_TKE_skewness_climatology-100m.png
export_fig -nocrop -png 20170101-20170520_arm-graciosa_TKE_skewness_climatology-100m.png
%%
monthnames = cellstr(datestr(datenum(2017,1,1):31:datenum(2017,12,31),'mmm'));
hf = figure; hf.Units = 'Normalized'; hf.Position = [0 0 1 1]; hf.Color = 'w';
for m = 1:12
    if strcmp(site,'arm-graciosa') && m < 6
        subplot(3,4,m);
        hold on
        
        %     p1 = patch([0:23 23:-1:0], [smooth(p25_latent.(['m' num2str(m)]),3,'rlowess'); flipud(smooth(p75_latent.(['m' num2str(m)]),3,'rlowess'))], rgb('DarkGray'),'EdgeColor','none');
        %     alpha(0.4)
        ax1 = gca; % current axes
        ax1.YColor = rgb('Crimson');
        ax1.XTick = 0:3:24;
        ax1.XLim = [0 24];
        ax1.YLim = [-5 -2];
        ax1_pos = ax1.Position; % position of first axes
        ax1.Box = 'on';
        ax1.YTickLabels = [repmat('10^{',length(ax1.YTick(:)),1) num2str(ax1.YTick(:)) repmat('}',length(ax1.YTick(:)),1)];
        
        %     e1 = line((0:23),med_latent.(['m' num2str(m)]),'Parent',ax1);
        switch m
            case {1,2,3,4,5}
                e1 = line((0.5:23.5),tke_med_300m(:,m),'Parent',ax1);
            case {6,7,8,9,10}
                e1 = line((0.5:23.5),tke_med_300m(:,m),'Parent',ax1);
            case {11,12,13,14,15}
                e1 = line((0.5:23.5),tke_med_300m(:,m),'Parent',ax1);
        end
        
        set(e1,'LineStyle',':','Marker','.','MarkerSize',10,'Color',rgb('Crimson'));
        
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','right',...
            'Color','none');
        ax2.YLim = [-1 1];
        ax2.XTick = 0:3:24;
        ax2.XLim = [0 24];
        ax2.YColor = rgb('Black');
        %     p2 = patch([0:23 23:-1:0], [smooth(p25_sensible.(['m' num2str(m)]),3,'rlowess'); flipud(smooth(p75_sensible.(['m' num2str(m)]),3,'rlowess'))], rgb('DarkSalmon'),'EdgeColor','none');
        %     alpha(0.4)
        
        switch m
            case {1,2,3,4,5}
                e2 = line((0.5:23.5),skewn_med_300m(:,m),'Parent',ax2);
            case {6,7,8,9,10}
                e2 = line((0.5:23.5),skewn_med_300m(:,m),'Parent',ax2);
            case {11,12,13,14,15}
                e2 = line((0.5:23.5),skewn_med_300m(:,m),'Parent',ax2);
        end
        set(e2,'LineStyle',':','Marker','.','MarkerSize',10,'Color',rgb('Black'));
        
        grid on;
        switch m
            case {9,10,11,12}
                xlabel('Time of day (UTC)');
        end
        switch m
            case {1,5,9}
                ylabel(ax1,'TKE (m^{2} s^{-3})');
        end
        switch m
            case {4,8,12}
                ylabel(ax2,'Skewness (unitless)');
        end
        if strcmp(site,'arm-graciosa') && m == 5
            ylabel(ax2,'Skewness (unitless)');
        end
        hold off
        title(monthnames{m})
    end
end
[hl] = legend([e1 e2],{' TKE,','Skewness'},'Orientation','Horizontal','Box','On');
hl.Position = [.275 .96 .5 .02];

% export_fig -nocrop -png 20160907-20170922_arm-ascension_TKE_skewness_climatology-300m.png
export_fig -nocrop -png 20170101-20170520_arm-graciosa_TKE_skewness_climatology-300m.png

%%
monthnames = cellstr(datestr(datenum(2017,1,1):31:datenum(2017,12,31),'mmm'));
hf = figure; hf.Units = 'Normalized'; hf.Position = [0 0 1 1]; hf.Color = 'w';
for m = 1:12
    if strcmp(site,'arm-graciosa') && m < 6
        subplot(3,4,m);
        hold on
        
        %     p1 = patch([0:23 23:-1:0], [smooth(p25_latent.(['m' num2str(m)]),3,'rlowess'); flipud(smooth(p75_latent.(['m' num2str(m)]),3,'rlowess'))], rgb('DarkGray'),'EdgeColor','none');
        %     alpha(0.4)
        ax1 = gca; % current axes
        ax1.YColor = rgb('Crimson');
        ax1.XTick = 0:3:24;
        ax1.XLim = [0 24];
        ax1.YLim = [-5 -2];
        ax1_pos = ax1.Position; % position of first axes
        ax1.Box = 'on';
        ax1.YTickLabels = [repmat('10^{',length(ax1.YTick(:)),1) num2str(ax1.YTick(:)) repmat('}',length(ax1.YTick(:)),1)];
        %     e1 = line((0:23),med_latent.(['m' num2str(m)]),'Parent',ax1);
        switch m
            case {1,2,3,4,5}
                e1 = line((0.5:23.5),tke_med_800m(:,m),'Parent',ax1);
            case {6,7,8,9,10}
                e1 = line((0.5:23.5),tke_med_800m(:,m),'Parent',ax1);
            case {11,12,13,14,15}
                e1 = line((0.5:23.5),tke_med_800m(:,m),'Parent',ax1);
        end
        
        set(e1,'LineStyle',':','Marker','.','MarkerSize',10,'Color',rgb('Crimson'));
        
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','bottom',...
            'YAxisLocation','right',...
            'Color','none');
        ax2.YLim = [-1 1];
        ax2.XTick = 0:3:24;
        ax2.XLim = [0 24];
        ax2.YColor = rgb('Black');
        %     p2 = patch([0:23 23:-1:0], [smooth(p25_sensible.(['m' num2str(m)]),3,'rlowess'); flipud(smooth(p75_sensible.(['m' num2str(m)]),3,'rlowess'))], rgb('DarkSalmon'),'EdgeColor','none');
        %     alpha(0.4)
        
        switch m
            case {1,2,3,4,5}
                e2 = line((0.5:23.5),skewn_med_800m(:,m),'Parent',ax2);
            case {6,7,8,9,10}
                e2 = line((0.5:23.5),skewn_med_800m(:,m),'Parent',ax2);
            case {11,12,13,14,15}
                e2 = line((0.5:23.5),skewn_med_800m(:,m),'Parent',ax2);
        end
        set(e2,'LineStyle',':','Marker','.','MarkerSize',10,'Color',rgb('Black'));
        
        grid on;
        switch m
            case {9,10,11,12}
                xlabel('Time of day (UTC)');
        end
        switch m
            case {1,5,9}
                ylabel(ax1,'TKE (m^{2} s^{-3})');
        end
        switch m
            case {4,8,12}
                ylabel(ax2,'Skewness (unitless)');
        end
        if strcmp(site,'arm-graciosa') && m == 5
            ylabel(ax2,'Skewness (unitless)');
        end
        hold off
        title(monthnames{m})
    end
    [hl] = legend([e1 e2],{' TKE,','Skewness'},'Orientation','Horizontal','Box','On');
    hl.Position = [.275 .96 .5 .02];
end
% export_fig -nocrop -png 20160907-20170922_arm-ascension_TKE_skewness_climatology-800m.png
export_fig -nocrop -png 20170101-20170520_arm-graciosa_TKE_skewness_climatology-800m.png

%%
% sp7 = subplot(3,3,7);
% % eb = errorbar(0.5:23.5,tke_med_100m(:,3),tke_p25_100m(:,3),tke_p75_100m(:,3));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% plot(.5:23.5,tke_med_100m(:,3),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_100m(:,3),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[2 1.2 8.5 3.8]); yticks = get(gca,'YTick');
% yticklabels = [repmat('10^{',length(yticks),1) num2str(yticks(:)) repmat('}',length(yticks),1)];
% set(gca,'YTickLabels',yticklabels,'Box','on'); 
% xlabel('Time UTC'); ylabel('TKE (m^{2} s^{-3})');
% text(.5,-2.2,'100 m')
% 
% sp8 = subplot(3,3,8);
% % eb = errorbar(0.5:23.5,tke_med_100m(:,4),tke_p25_100m(:,4),tke_p75_100m(:,4));
% hold on
% plot(.5:23.5,tke_med_100m(:,4),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_100m(:,4),'.','MarkerSize',10,'LineStyle',':','Color','r')
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[11.5 1.2 8.5 3.8],'YTickLabels', [],'Box','on'); 
% xlabel('Time UTC');
% text(.5,-2.2,'100 m')
% 
% sp9 = subplot(3,3,9);
% % eb = errorbar(0.5:23.5,tke_med_100m(:,5),tke_p25_100m(:,5),tke_p75_100m(:,5));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% plot(.5:23.5,tke_med_100m(:,5),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_100m(:,5),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[21 1.5 8.5 3.8],'YTickLabels', []); 
% xlabel('Time UTC');
% 
% sp4 = subplot(3,3,4);
% % eb = errorbar(0.5:23.5,tke_med_300m(:,3),tke_p25_300m(:,3),tke_p75_300m(:,3));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% plot(.5:23.5,tke_med_300m(:,3),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_300m(:,3),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[2 5.7 8.5 3.8],'XTickLabels',[]); yticks = get(gca,'YTick');
% yticklabels = [repmat('10^{',length(yticks),1) num2str(yticks(:)) repmat('}',length(yticks),1)];
% set(gca,'YTickLabels',yticklabels); 
% ylabel('TKE (m^{2} s^{-3})');
% text(.5,-2.2,'300 m')
% 
% sp5 = subplot(3,3,5);
% % eb = errorbar(0.5:23.5,tke_med_300m(:,4),tke_p25_300m(:,4),tke_p75_300m(:,4));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% plot(.5:23.5,tke_med_300m(:,4),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_300m(:,4),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[11.5 5.7 8.5 3.8],'YTickLabels', [],'XTickLabels',[]); 
% text(.5,-2.2,'300 m')
% 
% sp6 = subplot(3,3,6);
% % eb = errorbar(0.5:23.5,tke_med_300m(:,5),tke_p25_300m(:,5),tke_p75_300m(:,5));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% plot(.5:23.5,tke_med_300m(:,5),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_300m(:,5),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[21 6 8.5 3.8],'YTickLabels', [],'XTickLabels',[]); 
% 
% sp1 = subplot(3,3,1);
% % eb = errorbar(0.5:23.5,tke_med_800m(:,3),tke_p25_800m(:,3),tke_p75_800m(:,3));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r'; 
% plot(.5:23.5,tke_med_800m(:,3),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_800m(:,3),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[2 10.2 8.5 3.8],'XTickLabels',[]); yticks = get(gca,'YTick');
% yticklabels = [repmat('10^{',length(yticks),1) num2str(yticks(:)) repmat('}',length(yticks),1)];
% set(gca,'YTickLabels',yticklabels); 
% ylabel('TKE (m^{2} s^{-3})');
% text(.5,-2.2,'800 m'); title('March')
% 
% sp2 = subplot(3,3,2);
% % eb = errorbar(0.5:23.5,tke_med_800m(:,4),tke_p25_800m(:,4),tke_p75_800m(:,4));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% plot(.5:23.5,tke_med_800m(:,4),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_800m(:,4),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[11.5 10.2 8.5 3.8],'YTickLabels', [],'XTickLabels',[]); 
% text(.5,-2.2,'800 m'); title('April')
% 
% sp3 = subplot(3,3,3);
% % eb = errorbar(0.5:23.5,tke_med_800m(:,5),tke_p25_800m(:,5),tke_p75_800m(:,5));
% % eb.Marker = '.'; eb.MarkerSize = 10; eb.LineStyle = ':'; eb.Color = 'r';
% plot(.5:23.5,tke_med_800m(:,5),'LineStyle',':','Color','r')
% plot(.5:23.5,tke_med_800m(:,5),'.','MarkerSize',10,'LineStyle',':','Color','r')
% axis([0 24 -5 -2]); grid on;
% set(gca,'YTick',-6:-1,'XTick',0:4:24,'Units','centimeters','Position',[21 10.5 8.5 3.8],'YTickLabels', [],'XTickLabels',[]); 
% title('May')
% % export_fig -nocrop -png -painters 20170301-20170430_arm-graciosa_TKE-100-300-800m.png
% %%
% 
end
    


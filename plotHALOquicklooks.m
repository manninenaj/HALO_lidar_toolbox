function plotHALOquicklooks(site,DATES,processlev,measmode,typeof,varargin)
% Check inputs

p.ylabel = '';
p.masking = 0;
p.ylim = [0 99999];
p.ystep = 1000;
p.azilim = [0 360];
p.azistep = 60;

if ~isempty(varargin)
    p = parsePropertyValuePairs(p, varargin);
    if p.ylim(1)~=0
        p.ymin = p.ylim(1)/1000;
    end
    if p.ylim(2)~=99999
        p.ymax = p.ylim(2)/1000;
    end
    p.ystep = p.ystep/1000;
    if p.azilim(1)~=0
        p.azimin = p.azilim(1);
    end
    if p.azilim(2)~=360
        p.azimax = p.azilim(2);
    end
end


if nargin < 4
    error(sprintf(['At least inputs ''site'', ''DATE'', ''processlev'', and ''measmode'''...
        ' are required for the products: \n''epsilon'', ''wstats'', ''wstats4precipfilter'', ''sigma2vad''',...
        '''windshear'', ''LLJ'', ''ABLclassification'', ''cloud''']))
end
if nargin == 4 && (strcmp(processlev,'product') && any(strcmp(measmode,{'epsilon',...
        'wstats','wstats4precipfilter','sigma2vad','windshear','LLJ','ABLclassification','cloud'})) || ...
        strcmp(processlev,'background'))
    if ~ischar(site)
        error('The 1st input ''site'' must be a string.')
    end
    if length(DATES)>2
        error('''DATES'' can have max. length of 2.')
    elseif length(DATES)==1
        DATEstart = DATES; DATEend = DATES;
    elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
            length(num2str(DATES(2)))~=8)
        error(['The value(s) in the second input ''DATES'' must be' ...
            ' numerical date(s) in YYYYMMDD format.'])
    else
        DATEstart = DATES(1); DATEend = DATES(2);
    end
    if ~ischar(processlev) || ~any(strcmp(processlev,{'original','corrected','calibrated','background','product'}))
        error(['The 3rd input ''processlev'' must be a string and can be:'...
            ' ''original'', ''corrected'', ''calibrated'', ''background'', or ''product''.'])
    end
    if ~ischar(measmode) || ~any(strcmp(measmode,{'stare','vad','dbs','rhi','custom','co','cross','sector','windvad','winddbs',...
            'txt','wstats','wstats4precipfilter','epsilon','sigma2vad','windshear','LLJ','ABLclassification','cloud'}))
        error(sprintf(['The 4th input ''measmode'' must be a string and can be:\n'...
            '''stare'',''vad'',''dbs'',''rhi'',''co'',''cross'',''sector'',''custom'',''windvad'',''winddbs'',''txt'',''wstats''\n'...
            '''wstats4precipfilter'',''epsilon'',''sigma2vad'',''windshear'',''LLJ'',''ABLclassification'',''cloud''.']))
    end
end
if nargin < 5 && (~strcmp(processlev,'product') && ~any(strcmp(measmode,{'epsilon',...
        'wstats','wstats4precipfilter','sigma2vad','windshear','LLJ','ABLclassification','cloud'})) && ...
        ~strcmp(processlev,'background'))
    error(sprintf(['Inputs ''site'', ''DATE'', ''processlev'', ''measmode'', and ''typeof'''...
        ' are required for ANY OTHER products than: \n''epsilon'', ''wstats'', ''wstats4precipfilter'', ''sigma2vad'','...
        ' ''windshear'', ''LLJ'', ''ABLclassification'', ''cloud''']))
end
if nargin == 4 && (~strcmp(processlev,'product') && ~any(strcmp(measmode,{'epsilon',...
        'wstats','wstats4precipfilter','sigma2vad','windshear','LLJ','ABLclassification','cloud'})) && ...
        ~strcmp(processlev,'background'))
    error(sprintf(['Inputs ''site'', ''DATE'', ''processlev'', ''measmode'', and ''typeof'''...
        ' are required for ANY OTHER products than: \n''epsilon'', ''wstats'', ''wstats4precipfilter'', ''sigma2vad'','...
        ' ''windshear'', ''LLJ'', ''ABLclassification'', ''cloud''']))
else
    if ~ischar(site)
        error('The 1st input ''site'' must be a string.')
    end
    if length(DATES)>2
        error('''DATES'' can have max. length of 2.')
    elseif length(DATES)==1
        DATEstart = DATES; DATEend = DATES;
    elseif ~isnumeric(DATES) || (length(num2str(DATES(1)))~=8 && ...
            length(num2str(DATES(2)))~=8)
        error(['The value(s) in the second input ''DATES'' must be' ...
            ' numerical date(s) in YYYYMMDD format.'])
    else
        DATEstart = DATES(1); DATEend = DATES(2);
    end
    if ~ischar(processlev) || ~any(strcmp(processlev,{'original','corrected','calibrated','background','product'}))
        error(['The 3rd input ''processlev'' must be a string and can be:'...
            ' ''original'', ''corrected'', ''calibrated'', ''background'', or ''product''.'])
    end
    if ~ischar(measmode) || ~any(strcmp(measmode,{'stare','vad','dbs','rhi','co','cros','sector','custom','windvad','winddbs',...
            'txt','wstats','wstats4precipfilter','epsilon','sigma2vad','windshear','LLJ','ABLclassification','cloud'}))
        error(sprintf(['The 4th input ''measmode'' must be a string and can be:\n'...
            '''stare'',''vad'',''rhi'',''dbs'',''co'',''sector'',''cross'',''custom'',''windvad'',''winddbs'',''txt'',''wstats''\n'...
            '''wstats4precipfilter'',''epsilon'',''sigma2vad'',''windshear'',''LLJ'',''ABLclassification'',''cloud''.']))
    end
end

% Use datenum to accommodate leap years etc.
for DATEi = datenum(num2str(DATEstart),'yyyymmdd'):...
        datenum(num2str(DATEend),'yyyymmdd')
    
    % Convert date into required formats
    thedate = datestr(DATEi,'yyyymmdd');
    DATE = str2double(thedate);
    if exist('typeof','var') == 1
        [dirto,files] = getHALOfileList(site,DATE,processlev,measmode,typeof);
    else
        [dirto,files] = getHALOfileList(site,DATE,processlev,measmode);
    end
    if isempty(files), continue; end
    
    
    % Get default and site/unit/period specific parameters
    C = getconfig(site,DATE);
    
    switch processlev
        case 'calibrated'
            switch measmode
                case 'stare'
                    switch typeof
                        case {'co','cross','co12'}
                            data = load_nc_struct(fullfile([dirto files{1}]));
                            
                            if p.ymax == 99999
                                p.ymax = ceil(data.range(end)/1000);
                            end
                            tmax = p.ymax+p.ymax*.075;
                            hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                            hf.Color = 'white'; hf.Visible = 'off';
                            sp1 = subplot(321);
                            s0 = data.signal0;
                            imagesc(data.time,data.range/1000,transpose(s0)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[1.2 7.3 11 2.2],'Color',[.75 .75 .75]);
                            caxis([.99 1.02]); colormap(sp1,chilljet);
                            cb = colorbar; cb.Label.String = '(SNR+1)'; text(0,tmax,'uncorrected signal');
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 10.5; cb.Ticks = .99:.01:1.02; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)')
                            
                            sp2 = subplot(322);
                            s1 = data.signal;
                            imagesc(data.time,data.range/1000,transpose(s1)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[13.5 7.3 11 2.2],'Color',[.75 .75 .75]);
                            caxis([.99 1.02]); colormap(sp2,chilljet);
                            cb = colorbar; cb.Label.String = '(SNR+1)'; text(0,tmax,'corrected signal')
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 22.8;  cb.Ticks = .99:.01:1.02; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)');
                            
                            sp3 = subplot(323);
                            imagesc(data.time,data.range/1000,transpose(real(log10(data.beta_raw)))); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[1.2 4.2 11 2.2],'Color',[.75 .75 .75]);
                            caxis([-7 -4]); colormap(sp3,chilljet);
                            cb = colorbar; cb.Label.String = 'm-1 sr-1'; text(0,tmax,'beta');
                            cb.Ticks = -8:-3; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                                num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 10.5; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)')
                            
                            sp4 = subplot(324);
                            imagesc(data.time,data.range/1000,transpose(data.beta_error)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[13.5 4.2 11 2.2],'Color',[.75 .75 .75]);
                            caxis([0 1]); colormap(sp4,chilljet);
                            cb = colorbar; cb.Label.String = 'Fraction'; text(0,tmax,'signal error')
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 22.8; cb.Ticks = 0:.2:1; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)');
                            
                            sp5 = subplot(325);
                            imagesc(data.time,data.range/1000,transpose(data.v_raw)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[1.2 1.1 11 2.2],'Color',[.75 .75 .75]);
                            caxis([-2 2]); colormap(sp5,cmocean('balance'));
                            cb = colorbar; cb.Label.String = 'm s-1'; text(0,tmax,'vertical velocity')
                            ax1 = get(gca,'Position'); cb.Ticks = -2:.5:2; cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 10.5; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)'); xlabel('Time UTC (hours)')
                            
                            sp6 = subplot(326);
                            imagesc(data.time,data.range/1000,transpose(data.v_error)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[13.5 1.1 11 2.2],'Color',[.75 .75 .75]);
                            caxis([0 .5]); colormap(sp6,chilljet);
                            cb = colorbar; cb.Label.String = 'm s-1'; cb.Ticks = 0:.1:.5; text(0,tmax,'vertical velocity error')
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            xlabel('Time UTC (hours)');ylabel('Range (km)');
                            
                            fname = fullfile([dirto strrep(files{1},'.nc','.png')]);
                            fprintf('Writing %s\n',fname)
                            export_fig('-png','-m2',fname)
                            
                            close(hf)
                        otherwise
                            continue
                    end
                case {'vad','sector'}
                    
                    for i = 1:length(files)
                        %if i~=1, return; end
                        data = load_nc_struct(fullfile([dirto files{i}]));
                        if p.masking~=0
                            smask = logical(zeros(size(data.signal)));
                            %for ia = 1:length(data.azimuth)
                            %   ibegin = find(data.signal(ia,5:end)>1.2,1,'first');
                            %   smask(ia,ibegin+4:end) = true;
                            %end
                            smask(data.signal<p.masking) = true;
                            smask(:,1:4) = true;
                            data.beta_raw(smask) = nan;
                            data.v_raw(smask) = nan;
                            data.signal(smask) = nan;
                            data.beta_error(smask) = nan;
                            data.v_error(smask) = nan;
                        end
                        
                        if p.ymax~=99999 % default
                            cond = data.range/1000>p.ymax;
                            
                            data.signal(:,cond) = [];
                            data.beta_raw(:,cond) = [];
                            data.v_raw(:,cond) = [];
                            data.beta_error(:,cond) = [];
                            data.v_error(:,cond) = [];
                            data.range(cond) = [];
                        end
                        r = transpose(data.range(:)/1000);
                        a = transpose(data.azimuth(:));
                        a(a>360) = a(a>360)-360;
                        s = transpose(10*real(log10(data.signal-1)));
                        s_e = transpose(data.beta_error * 100);
                        b = transpose(real(log10(data.beta_raw)));
                        b_e = transpose(data.beta_error * 100);
                        v = transpose(data.v_raw);
                        v_e = transpose(data.v_error);
                        
                        
                        hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 30 12];
                        hf.Color = 'white'; hf.Visible = 'off';
                        sp1 = subplot(231);;
                        [~,c]= polarPcolor(r,a,s,'rStep',p.ystep,'thetaStep',p.azistep,'labelR',p.ylabel,'thetaMin',p.azimin,'thetaMax',p.azimax);
                        ylabel(c,' signal intensity (dB)');
                        set(gcf,'color','w')
                        colormap(sp1,flipud(cbrewer('div','Spectral',256)))
                        caxis([10*log10(p.masking-1) 0])
                        
                        sp2 = subplot(232);
                        [~,c]= polarPcolor(r,a,b,'rStep',p.ystep,'thetaStep',p.azistep,'labelR',p.ylabel,'thetaMin',p.azimin,'thetaMax',p.azimax);
                        ylabel(c,' log_{10} att. beta (m-1 sr-1)');
                        set(gcf,'color','w')
                        colormap(sp2,flipud(cbrewer('div','Spectral',256)))
                        caxis([-7 -4])
                        
                        sp3 = subplot(233);
                        [~,c]= polarPcolor(r,a,v,'rStep',p.ystep,'thetaStep',p.azistep,'labelR',p.ylabel,'thetaMin',p.azimin,'thetaMax',p.azimax);
                        ylabel(c,' radial velocity (m s-1)');
                        set(gcf,'color','w')
                        colormap(sp3,cmocean('balance'))
                        caxis([-20 20])
                        
                        sp4 = subplot(234);
                        [~,c]= polarPcolor(r,a,s_e,'rStep',p.ystep,'thetaStep',p.azistep,'labelR',p.ylabel,'thetaMin',p.azimin,'thetaMax',p.azimax);
                        ylabel(c,' signal error (%)');
                        set(gcf,'color','w')
                        colormap(sp4,flipud(cbrewer('div','Spectral',256)))
                        caxis([0 100])
                        
                        sp5 = subplot(235);
                        [~,c]= polarPcolor(r,a,b_e,'rStep',p.ystep,'thetaStep',p.azistep,'labelR',p.ylabel,'thetaMin',p.azimin,'thetaMax',p.azimax);
                        ylabel(c,' att. beta error (%)');
                        set(gcf,'color','w')
                        colormap(sp5,flipud(cbrewer('div','Spectral',256)))
                        caxis([0 100])
                        
                        sp6 = subplot(236);
                        [~,c]= polarPcolor(r,a,v_e,'rStep',p.ystep,'thetaStep',p.azistep,'labelR',p.ylabel,'thetaMin',p.azimin,'thetaMax',p.azimax);
                        ylabel(c,' radial velocity error (m s-1)');
                        set(gcf,'color','w')
                        colormap(sp6,flipud(cbrewer('div','Spectral',256)))
                        caxis([0 .5])
                        
                        set(findall(hf,'-property','FontSize'),'FontSize',8)
                        
                        fname = fullfile([dirto strrep(files{i},'.nc','.png')]);
                        fprintf('Writing \n    %s\n',fname)
                        export_fig('-png','-m2',fname)
                        
                        close(hf)
                        
                    end
                    
                otherwise
                    continue
            end
            
        case 'original'
            switch measmode
                case 'stare'
                    switch typeof
                        case 'co'
                            data = load_nc_struct(fullfile([dirto files{1}]));
                            
                            tmax = p.ymax+p.ymax*.075;
                            hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                            hf.Color = 'white'; hf.Visible = 'off';
                            sp1 = subplot(321);
                            s0 = data.signal;
                            imagesc(data.time,data.range/1000,transpose(s0)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[1 7.3 11 2.2],'Color',[.75 .75 .75]);
                            caxis([.99 1.01]); colormap(sp1,flipud(cbrewer('div','Spectral',256)));
                            cb = colorbar; cb.Label.String = '(SNR+1)'; text(0,tmax,'uncorrected signal');
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)')
                            
                            sp2 = subplot(322);
                            s1 = nan(size(data.signal));
                            imagesc(data.time,data.range/1000,transpose(s1)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[13.5 7.3 11 2.2],'Color',[.75 .75 .75]);
                            caxis([-30 0]); colormap(sp2,flipud(cbrewer('div','Spectral',256)));
                            cb = colorbar; cb.Label.String = 'dB'; text(0,tmax,'corrected signal')
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            X_out1 = [];
                            
                            sp3 = subplot(323);
                            imagesc(data.time,data.range/1000,transpose(real(log10(data.beta_raw)))); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[1 4.2 11 2.2],'Color',[.75 .75 .75]);
                            caxis([-7.1 -4]); colormap(sp3,flipud(cbrewer('div','Spectral',256)));
                            cb = colorbar; cb.Label.String = 'm-1 sr-1'; text(0,tmax,'beta');
                            cb.Ticks = -8:-3; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                                num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)')
                            X_out1 = [];
                            
                            sp4 = subplot(324);
                            imagesc(data.time,data.range/1000,transpose(data.beta_error)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[13.5 4.2 11 2.2],'Color',[.75 .75 .75]);
                            caxis([0 1]); colormap(sp4,flipud(cbrewer('div','Spectral',256)));
                            cb = colorbar; cb.Label.String = 'Fraction'; text(0,tmax,'beta error')
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            X_out1 = [];
                            
                            sp5 = subplot(325);
                            imagesc(data.time,data.range/1000,transpose(data.v_raw)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[1 1.1 11 2.2],'Color',[.75 .75 .75]);
                            caxis([-1.5 1.5]); colormap(sp5,cmocean('balance'));
                            cb = colorbar; cb.Label.String = 'm s-1'; text(0,tmax,'vertical velocity')
                            ax1 = get(gca,'Position'); cb.Ticks = -1.5:.5:1.5; cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            ylabel('Range (km)'); xlabel('Time UTC')
                            X_out1 = [];
                            
                            sp6 = subplot(326);
                            imagesc(data.time,data.range/1000,transpose(data.v_error)); axis([0 24 0 p.ymax]); shading flat;
                            set(gca,'YDir','normal','Ytick',0:p.ystep:ceil(p.ymax),'XTick',0:3:24,'Units',...
                                'centimeters','Position',[13.5 1.1 11 2.2],'Color',[.75 .75 .75]);
                            caxis([0 .5]); colormap(sp6,flipud(cbrewer('div','Spectral',256)));
                            cb = colorbar; cb.Label.String = 'm s-1'; text(0,tmax,'vertical velocity error')
                            ax1 = get(gca,'Position'); cb.Units = 'centimeters'; cb.Position(3) = .25;
                            cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                            xlabel('Time UTC')
                            
                            fname = fullfile([dirto strrep(files{1},'.nc','.png')]);
                            fprintf('Writing %s\n',fname)
                            export_fig('-png','-m1',fname)
                            close(hf)
                            
                            
                        otherwise
                            continue
                    end
                case 'vad'
                    for i = 1:length(files)
                        data = load_nc_struct(fullfile([dirto files{i}]),{'range','azimuth','v_raw','beta_raw','signal'});
                        if p.masking
                            smask = logical(zeros(size(data.signal)));
                            for ia = 1:length(data.azimuth)
                                ibegin = find(data.signal(ia,5:end)>1.2,1,'first');
                                smask(ia,4+ibegin:end) = true;
                            end
                            smask(data.signal<1.01) = true;
                            smask(:,1:4) = true;
                            data.beta_raw(smask) = nan;
                            data.v_raw(smask) = nan;
                            data.signal(smask) = nan;
                        end
                        if p.ymax~=99999
                            cond = data.range>p.ymax;
                            data.signal(:,cond) = [];
                            data.beta_raw(:,cond) = [];
                            data.v_raw(:,cond) = [];
                            data.range(cond) = [];
                        end
                        r = transpose(data.range(:)/1000);
                        a = transpose(data.azimuth(:));
                        a(a>360) = a(a>360)-360;
                        s = transpose(real(log10(data.signal-1)));
                        s_e = transpose(data.beta_error * 100);
                        b = transpose(real(log10(data.beta_raw)));
                        b_e = transpose(data.beta_error * 100);
                        v = transpose(data.v_raw);
                        v_e = transpose(data.v_error);
                        
                        %circles = linspace(round(r(1),1),round(r(end),1),ncircles);
                        nspokes = 9;
                        %rticklabel = cellstr(num2str(circles(:)));
                        
                        hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 30 10];
                        hf.Color = 'white'; hf.Visible = 'off';
                        sp1 = subplot(231);
                        [~,c]= polarPcolor(r,a,s,'rStep',p.ystep,'Nspokes',nspokes);
                        ylabel(c,' signal intensity (dB)');
                        set(gcf,'color','w')
                        colormap(sp1,flipud(cbrewer('div','Spectral',256)))
                        caxis([-20 0])
                        
                        sp2 = subplot(232);
                        [~,c]= polarPcolor(r,a,b,'rStep',p.ystep,'Nspokes',nspokes);
                        ylabel(c,' log_{10} att. beta (m-1 sr-1)');
                        set(gcf,'color','w')
                        colormap(sp2,flipud(cbrewer('div','Spectral',256)))
                        caxis([-7 -4])
                        
                        sp3 = subplot(233);
                        [~,c]= polarPcolor(r,a,v,'rStep',p.ystep,'Nspokes',nspokes);
                        ylabel(c,' radial velocity (m s-1)');
                        set(gcf,'color','w')
                        colormap(sp3,cmocean('balance'))
                        caxis([-10 10])
                        
                        sp4 = subplot(234);
                        [~,c]= polarPcolor(r,a,s_e,'rStep',p.ystep,'Nspokes',nspokes);
                        ylabel(c,' signal error (%)');
                        set(gcf,'color','w')
                        colormap(sp4,flipud(cbrewer('div','Spectral',256)))
                        caxis([0 100])
                        
                        sp5 = subplot(235);
                        [~,c]= polarPcolor(r,a,b_e,'rStep',p.ystep,'Nspokes',nspokes);
                        ylabel(c,' att. beta error (%)');
                        set(gcf,'color','w')
                        colormap(sp5,flipud(cbrewer('div','Spectral',256)))
                        caxis([0 100])
                        
                        sp6 = subplot(236);
                        [~,c]= polarPcolor(r,a,v_e,'rStep',p.ystep,'Nspokes',nspokes);
                        ylabel(c,' radial velocity error (m s-1)');
                        set(gcf,'color','w')
                        colormap(sp6,flipud(cbrewer('div','Spectral',256)))
                        caxis([0 .5])
                        
                        set(findall(hf,'-property','FontSize'),'FontSize',8)
                        
                        fname = fullfile([dirto strrep(files{i},'.nc','.png')]);
                        fprintf('Writing %s\n',fname)
                        export_fig('-png','-m2',fname)
                        close(hf)
                        
                    end
                    
                otherwise
                    
                    continue
            end
        case 'product'
            switch measmode
                case 'wstats'
                    data = load_nc_struct(fullfile([dirto files{1}]));
                    
                    beta_mean = data.beta_mean_3min;
                    snr_mean = data.signal_mean_3min;
                    velo_mean = data.radial_velocity_mean_3min;
                    velo_var = data.radial_velocity_variance_3min;
                    if isfield(data,'radial_velocity_skewness_60min')
                        velo_skewn = data.radial_velocity_skewness_60min;
                        velo_kurto = data.radial_velocity_kurtosis_60min;
                        wstats_time = data.time_60min;
                    else
                        velo_skewn = data.radial_velocity_skewness_30min;
                        velo_kurto = data.radial_velocity_kurtosis_30min;
                        wstats_time = data.time_30min;
                    end
                    %% Create a cleaning filter based what field are available
                    %fnames = fieldnames(data);
                    %switch ~isempty(strmatch('signal_instrumental_precision_variance_',fnames))
                    %    case 1
                    %        condnan = ...
                    %            10*real(log10(data.signal_mean_3min-1)) < -23 | ...
                    %            isnan(data.signal_mean_3min) | ...
                    %            real(log10(data.signal_instrumental_precision_variance_3min)) > 0 | ...
                    %            data.nsamples_3min < round(max(data.nsamples_3min(:))*.66);
                    %    case 0
                    %        condnan = ...
                    %            10*real(log10(data.signal_mean_3min-1)) < -23 | ...
                    %            isnan(data.signal_mean_3min) | ...
                    %            data.nsamples_3min < round(max(data.nsamples_3min(:))*.66);
                    %end
                    %beta_mean(condnan) = nan;
                    %snr_var(condnan) = nan;
                    %velo_mean(condnan) = nan;
                    %velo_var(condnan) = nan;
                    %condnan = 10*real(log10(data.signal_mean_60min-1))<-20 | ...
                    %    isnan(data.signal_mean_60min) | ...
                    %    data.radial_velocity_mean_error_60min > .33;
                    %velo_skewn(condnan) = nan;
                    %velo_kurto(condnan) = nan;
                    
                    hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                    hf.Color = 'white'; hf.Visible = 'off';

                    sp1 = subplot(321);
                    pcolor(data.time_3min,data.height/1000,transpose(real(log10(beta_mean)))); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 7.3 11 2.2]);
                    caxis([-7 -4]); colormap(sp1,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'beta mean');
                    cb = colorbar; cb.Label.String = 'm-1 sr-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Ticks = -7:-4; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                        num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)');
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp2 = subplot(322);
                    pcolor(data.time_3min,data.height/1000,transpose(snr_mean)); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 7.3 11 2.2]);
                    caxis([0.995 1.01]); colormap(sp2,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'signal mean')
                    cb = colorbar; cb.Label.String = 'SNR+1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = .995:.005:1.01; ylabel('Height (km)');
                    set(gca,'Color',[.5 .5 .5])
                     
                    sp3 = subplot(323);
                    pcolor(data.time_3min,data.height/1000,transpose(velo_mean)); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 4.2 11 2.2]);
                    caxis([-3 3]); colormap(sp3,cmocean('balance')); text(0,3.35,'velocity mean');
                    cb = colorbar; cb.Label.String = 'm s-1'; cb.Ticks = -3:1:3; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)')
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp4 = subplot(324);
                    pcolor(data.time_3min,data.height/1000,transpose(real(log10(velo_var)))); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 4.2 11 2.2]);
                    caxis([-3 1]); colormap(sp4,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'velocity variance')
                    cb = colorbar; cb.Label.String = 'm2 s-2'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = -3:1; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                        num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    ylabel('Height (km)')
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp5 = subplot(325);
                    pcolor(wstats_time,data.height/1000,transpose(velo_skewn)); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 1.1 11 2.2]);
                    caxis([-2 2]); colormap(sp5,cmocean('balance')); text(0,3.35,'velocity skewness')
                    cb = colorbar; cb.Label.String = '-'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)'); xlabel('Time UTC'); cb.Ticks = -2:1:2;
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp6 = subplot(326);
                    pcolor(wstats_time,data.height/1000,transpose(velo_kurto)); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 1.1 11 2.2]);
                    caxis([-4 6]); colormap(sp6,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'velocity kurtosis')
                    cb = colorbar; cb.Label.String = '-'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25;   cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)'); xlabel('Time UTC')
                    set(gca,'Color',[.5 .5 .5])
                    
                    if exist('typeof','var') == 1
                        [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode,typeof);
                        export_fig('-png',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                            '-%s-%s.png'], dir_out,num2str(DATE),site,measmode,typeof))
                    else
                        [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode);
                        export_fig('-png','-m2',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                            '-%s.png'], dir_out,num2str(DATE),site,measmode))
                    end
                    close(hf)
                case 'epsilon'
                    data = load_nc_struct(fullfile([dirto files{1}]));
                    
                    [dirsnr,filessnr] = getHALOfileList(site,DATE,'product' ,'wstats');
                    wstats = load_nc_struct(fullfile([dirsnr filessnr{1}]));
                    
                    epsilon = data.epsilon_3min;
                    epsilon_error = data.epsilon_error_3min;
                    L = data.L_3min;
                    L1 = data.L1_3min;
                    velovar_error = wstats.radial_velocity_simple_variance_error_3min;
                    velo_var = wstats.radial_velocity_simple_variance_3min - wstats.radial_velocity_instrumental_precision_variance_3min;
                    velo_var(velo_var<0)=nan;
                    
                    cond = [];
                    %10*real(log10(wstats.signal_mean_3min-1))<-32 | ...
                    %isnan(wstats.signal_mean_3min);% | ...
                    %                        wstats.radial_velocity_mean_error_3min > .33;
                    epsilon(cond) = nan;
                    epsilon_error(cond) = nan;
                    L(cond) = nan;
                    L1(cond) = nan;
                    velovar_error(cond) = nan;
                    velo_var(cond) = nan;
                    
                    hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                    hf.Color = 'white'; hf.Visible = 'off';
                    sp1 = subplot(321);
                    pcolor(data.time_3min,data.height/1000,real(log10(epsilon))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 7.3 11 2.2],'Color',rgb('DarkGray'));
                    caxis([-7 -1]); colormap(sp1,morgenstemning); text(0,3.35,'epsilon');
                    cb = colorbar; cb.Label.String = 'm2 s-3'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Ticks = -7:-1; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                        num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)')
                    sp2 = subplot(322);
                    pcolor(data.time_3min,data.height/1000,epsilon_error'*100); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 7.3 11 2.2],'Color',rgb('DarkGray'));
                    caxis([0 300]); colormap(sp2,cmap_darkviolet_to_brickred); text(0,3.35,'epsilon error')
                    cb = colorbar; cb.Label.String = '%'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.7; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)')
                    sp3 = subplot(323);
                    pcolor(data.time_3min,data.height/1000,real(log10(L))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 4.2 11 2.2],'Color',rgb('DarkGray'));
                    caxis([2 4]); colormap(sp3,cmap_darkviolet_to_brickred); text(0,3.35,'L');
                    cb = colorbar; cb.Label.String = 'm'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = 2:4; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                        num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    ylabel('Height (km)')
                    sp4 = subplot(324);
                    pcolor(data.time_3min,data.height/1000,real(log10(L1))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 4.2 11 2.2],'Color',rgb('DarkGray'));
                    caxis([1 3]); colormap(sp4,cmap_darkviolet_to_brickred); text(0,3.35,'L1')
                    cb = colorbar; cb.Label.String = 'm'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.7; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = 1:3; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                        num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    ylabel('Height (km)')
                    sp5 = subplot(325);
                    pcolor(data.time_3min,data.height/1000,real(log10(velo_var))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 1.1 11 2.2],'Color',rgb('DarkGray'));
                    caxis([-4 1]); colormap(sp5,cmap_darkviolet_to_brickred); text(0,3.35,'velocity variance')
                    cb = colorbar; cb.Label.String = 'm2 s-2'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = -4:1; cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                        num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    ylabel('Height (km)')
                    sp6 = subplot(326);
                    pcolor(wstats.time_3min,wstats.height/1000,real(log10(velovar_error))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 1.1 11 2.2],'Color',rgb('DarkGray'));
                    caxis([-3 0]); colormap(sp6,cmap_darkviolet_to_brickred); text(0,3.35,'velocity variance error')
                    cb = colorbar; cb.Label.String = 'm2 s-2'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25;   cb.Position(1) = 22.7; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)'); xlabel('Time UTC')
                    
                    [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode);
                    export_fig('-png','-m2',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                        '-%s.png'], dir_out,num2str(DATE),site,measmode))
                    close(hf)
                case 'windshear'
                    data = load_nc_struct(fullfile([dirto files{1}]));
                    
                    
                    windhsear = data.vector_wind_shear_3min;
                    windhsear_e = data.vector_wind_shear_error_3min;
                    
                    hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                    hf.Color = 'white'; hf.Visible = 'off';
                    sp1 = subplot(321);
                    pcolor(data.time_3min,data.height/1000,((windhsear))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 7.3 11 2.2]);
                    caxis([0 0.06]); colormap(sp1,cmocean('thermal')); text(0,3.35,'vector wind shear');
                    cb = colorbar; cb.Label.String = 's-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Ticks = 0:0.01:0.06;
                    %                     cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                    %                         num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    
                    sp1 = subplot(322);
                    pcolor(data.time_3min,data.height/1000,((windhsear_e))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 4.2 11 2.2]);
                    caxis([0 0.01]); colormap(sp1,cmocean('thermal')); text(0,3.35,'vector wind shear error');
                    cb = colorbar; cb.Label.String = 's-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Ticks = 0:0.01:0.06;
                    %                     cb.TickLabels = [repmat('10^{',length(cb.Ticks(:)),1), ...
                    %                         num2str(cb.Ticks(:)) repmat('}',length(cb.Ticks(:)),1)];
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    
                    
                    [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode);
                    export_fig('-png','-m2',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                        '-%s.png'], dir_out,num2str(DATE),site,measmode))
                    close(hf)
                case 'ABLclassification'
                    [dir_bl,files_bl] = getHALOfileList(site,DATE,'product' ,'ABLclassification');
                    [bl,blatt] = load_nc_struct(fullfile([dir_bl files_bl{1}]));
                    
                    TKEconnected = double(bl.turbulence_coupling_3min);
                    TKEconnected(TKEconnected==0) = nan;
                    BLclass = double(bl.bl_classification_3min);
                    BLclass(BLclass==0) = nan;
                    
                    cmap_tkecw = [blatt.turbulence_coupling_3min.legend_key_red;...
                        blatt.turbulence_coupling_3min.legend_key_green;...
                        blatt.turbulence_coupling_3min.legend_key_blue]';
                    
                    cmap_blc = [blatt.bl_classification_3min.legend_key_red;...
                        blatt.bl_classification_3min.legend_key_green;...
                        blatt.bl_classification_3min.legend_key_blue]';
                    
                    hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                    hf.Color = 'white'; hf.Visible = 'off';
                    
                    sp1 = subplot(321);
                    pcolor(bl.time_3min,bl.height/1000,TKEconnected'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 7.3 11 2.2]);
                    caxis([0 7]); colormap(sp1,cmap_tkecw); text(0,3.35,'Turbulence coupling');
                    cb = colorbar; cb.Ticks = .5:6.5;  cb.Units = 'centimeters';
                    cb.TickLabels = sprintf(blatt.turbulence_coupling_3min.definition); cb.FontSize = 8; pause(.5); ax1 = get(gca,'Position');
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)');
                    
                    sp3 = subplot(323);
                    pcolor(bl.time_3min,bl.height/1000,BLclass'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 4.2 11 2.2]);
                    caxis([0 10]); colormap(sp3,cmap_blc); text(0,3.35,'Boundary layer classification')
                    cb = colorbar; cb.Ticks = .5:9.5;  cb.Units = 'centimeters';
                    cb.TickLabels = sprintf(blatt.bl_classification_3min.definition); cb.FontSize = 8; pause(.5); ax1 = get(gca,'Position');
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)'); xlabel('Time UTC'); pause(.1)
                    
                    [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode);
                    export_fig('-png','-m2',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                        '-%s.png'], dir_out,num2str(DATE),site,measmode))
                    close(hf)
                    
                case 'windvad'
                    [data] = load_nc_struct(fullfile([dirto files{1}]));
                    
                    hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                    hf.Color = 'white'; hf.Visible = 'off';
                    
                    snr = transpose(data.mean_snr);
                    ws = transpose(data.wind_speed);
                    wd = transpose(data.wind_direction);
                    ws_e = transpose(data.wind_speed_error);
                    wd_e = transpose(data.wind_direction_error);
                    w = transpose(data.w);
                                       
                    if ~isfield(data,'height_agl')
                        height = data.height/1000;
                    else
                        height = data.height_agl/1000;
                    end
                    hlabel = 'Height agl (km)';
                    
                    sp1 = subplot(321);
                    pcolor(data.time,height,ws); axis([0 24 0 p.ymax]); shading flat
                    set(gca,'XTick',0:3:24,'Units','centimeters','Position',[1 7.3 11 2.2],'Color',rgb('DarkGray'),'YTick',0:p.ystep:ceil(p.ymax));
                    caxis([0 20]); colormap(sp1,cmocean('thermal')); text(0,p.ymax+p.ymax*.1,'Wind speed');
                    cb = colorbar; cb.Label.String = 'm s-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Ticks = 0:5:20; cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel(hlabel);
                    
                    sp2 = subplot(322);
                    pcolor(data.time,height,ws_e); axis([0 24 0 p.ymax]); shading flat
                    set(gca,'XTick',0:3:24,'Units','centimeters','Position',[13.5 7.3 11 2.2],'Color',rgb('DarkGray'),'YTick',0:p.ystep:ceil(p.ymax));
                    caxis([0 3]); colormap(sp2,cmocean('thermal')); text(0,p.ymax+p.ymax*.1,'Wind speed error')
                    cb = colorbar; cb.Ticks = 0:.5:10; cb.Label.String = 'm s-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel(hlabel);
                    
                    sp3 = subplot(323);
                    pcolor(data.time,height,wd); axis([0 24 0 p.ymax]); shading flat
                    set(gca,'XTick',0:3:24,'Units','centimeters','Position',[1 4.2 11 2.2],'Color',rgb('DarkGray'),'YTick',0:p.ystep:ceil(p.ymax));
                    caxis([0 360]); colormap(sp3,colorcet('C8')); text(0,p.ymax+p.ymax*.1,'Wind direction');
                    cb = colorbar; cb.Label.String = 'degrees'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = 0:90:360;
                    ylabel(hlabel);
                    
                    sp4 = subplot(324);
                    pcolor(data.time,height,wd_e); axis([0 24 0 p.ymax]); shading flat
                    set(gca,'XTick',0:3:24,'Units','centimeters','Position',[13.5 4.2 11 2.2],'Color',rgb('DarkGray'),'YTick',0:p.ystep:ceil(p.ymax));
                    caxis([0 2]); colormap(sp4,cmocean('thermal')); text(0,p.ymax+p.ymax*.1,'Wind direction error')
                    cb = colorbar; cb.Label.String = 'degrees'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = 0:.5:2;
                    ylabel(hlabel);
                    
                    sp5 = subplot(325);
                    pcolor(data.time,height,w); axis([0 24 0 p.ymax]); shading flat
                    set(gca,'XTick',0:3:24,'Units','centimeters','Position',[1 1.1 11 2.2],'Color',rgb('DarkGray'),'YTick',0:p.ystep:ceil(p.ymax));
                    caxis([-3 3]); colormap(sp5,cmocean('balance')); text(0,p.ymax+p.ymax*.1,'w wind component')
                    cb = colorbar; cb.Label.String = 'm s-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = -3:1:3;
                    ylabel(hlabel);
                    
                    sp6 = subplot(326);
                    pcolor(data.time,height,snr); axis([0 24 0 p.ymax]); shading flat
                    set(gca,'XTick',0:3:24,'Units','centimeters','Position',[13.5 1.1 11 2.2],'Color',rgb('DarkGray'),'YTick',0:p.ystep:ceil(p.ymax));
                    caxis([.995 1.015]); colormap(sp6,flipud(cbrewer('div','Spectral',256))); text(0,p.ymax+p.ymax*.1,'Mean signal')
                    cb = colorbar; cb.Label.String = 'SNR+1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = .995:.005:1.015;
                    ylabel(hlabel);
                    
                    [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode,typeof);
                    fname = strrep(files{1},'.nc','.png');
                    
                    fprintf('Writing %s\n',fullfile([dir_out fname]))
                    export_fig('-png','-m2',fullfile([dir_out fname]))
                    close(hf)

                    
                case 'winddbs'
                    [data] = load_nc_struct(fullfile([dirto files{1}]));
                    
                    hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                    hf.Color = 'white'; hf.Visible = 'off';
                    
                    snr = data.mean_snr;
                    ws = data.wind_speed;
                    wd = data.wind_direction;
                    %                     ws_e = data.wind_speed_error;
                    %                     wd_e = data.wind_direction_error;
                    w = data.w;
                    
                    sp1 = subplot(321);
                    pcolor(data.time,data.height/1000,ws'); axis([0 24 0 4]); shading flat
                    set(gca,'Ytick',0:4,'XTick',0:3:24,'Units','centimeters','Position',[1 7.3 11 2.2],'Color',rgb('DarkGray'));
                    caxis([0 20]); colormap(sp1,cmocean('thermal')); text(0,4.45,'Wind speed');
                    cb = colorbar; cb.Label.String = 'm s-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Ticks = 0:5:20; cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)');
                    
                    sp2 = subplot(322);
                    %                     pcolor(data.time,data.height/1000,ws_e');
                    axis([0 24 0 4]); shading flat
                    set(gca,'Ytick',0:4,'XTick',0:3:24,'Units','centimeters','Position',[13.5 7.3 11 2.2],'Color',rgb('DarkGray'));
                    caxis([0 3]); colormap(sp2,cmocean('thermal')); text(0,4.45,'Wind speed error')
                    cb = colorbar; cb.Ticks = 0:.5:10; cb.Label.String = 'm s-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)')
                    
                    sp3 = subplot(323);
                    pcolor(data.time,data.height/1000,wd'); axis([0 24 0 4]); shading flat
                    set(gca,'Ytick',0:4,'XTick',0:3:24,'Units','centimeters','Position',[1 4.2 11 2.2],'Color',rgb('DarkGray'));
                    caxis([0 360]); colormap(sp3,colorcet('C8')); text(0,4.45,'Wind direction');
                    cb = colorbar; cb.Label.String = 'degrees'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = 0:90:360;
                    ylabel('Height (km)')
                    
                    sp4 = subplot(324);
                    %                     pcolor(data.time,data.height/1000,wd_e');
                    axis([0 24 0 4]); shading flat
                    set(gca,'Ytick',0:4,'XTick',0:3:24,'Units','centimeters','Position',[13.5 4.2 11 2.2],'Color',rgb('DarkGray'));
                    caxis([0 2]); colormap(sp4,cmocean('thermal')); text(0,4.45,'Wind direction error')
                    cb = colorbar; cb.Label.String = 'degrees'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = 0:.5:2;
                    ylabel('Height (km)')
                    
                    sp5 = subplot(325);
                    pcolor(data.time,data.height/1000,w'); axis([0 24 0 4]); shading flat
                    set(gca,'Ytick',0:4,'XTick',0:3:24,'Units','centimeters','Position',[1 1.1 11 2.2],'Color',rgb('DarkGray'));
                    caxis([-3 3]); colormap(sp5,cmocean('balance')); text(0,4.45,'w wind component')
                    cb = colorbar; cb.Label.String = 'm s-1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 10.2; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = -3:1:3;
                    ylabel('Height (km)')
                    
                    sp6 = subplot(326);
                    pcolor(data.time,data.height/1000,snr'); axis([0 24 0 4]); shading flat
                    set(gca,'Ytick',0:4,'XTick',0:3:24,'Units','centimeters','Position',[13.5 1.1 11 2.2],'Color',rgb('DarkGray'));
                    caxis([.995 1.015]); colormap(sp6,cmap_darkviolet_to_brickred); text(0,4.45,'Mean signal')
                    cb = colorbar; cb.Label.String = 'SNR+1'; ax1 = get(gca,'Position'); cb.Units = 'centimeters';
                    cb.Position(3) = .25; cb.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb.Ticks = .995:.005:1.015;
                    ylabel('Height (km)')
                    
                    [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode,typeof);
                    export_fig('-png','-m2',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                        '-%s.png'], dir_out,num2str(DATE),site,measmode))
                    close(hf)
                    
                case 'wstats4precipfilter'
                    data = load_nc_struct(fullfile([dirto files{1}]));
                    
                    beta_mean = data.beta_mean_3min;
                    beta_var = data.beta_variance_3min;
                    snr_mean = data.signal_mean_3min;
                    snr_var = data.signal_variance_3min;
                    snr_instr_var = data.signal_instrumental_precision_variance_3min;
                    velo_mean = data.radial_velocity_mean_3min;
                    % Create a cleaning filter based what field are available
                    fnames = fieldnames(data);
                    switch ~isempty(strmatch('signal_instrumental_precision_variance_',fnames))
                        case 1
                            condnan = ...
                                10*real(log10(data.signal_mean_3min-1)) < -23 | ...
                                isnan(data.signal_mean_3min) | ...
                                data.nsamples_3min < round(max(data.nsamples_3min(:))*.66);% real(log10(data.signal_instrumental_precision_variance_3min)) > 0 | ...
                        case 0
                            condnan = ...
                                10*real(log10(data.signal_mean_3min-1)) < -23 | ...
                                isnan(data.signal_mean_3min) | ...
                                data.nsamples_3min < round(max(data.nsamples_3min(:))*.66);
                    end
                    condnan(:,1:3) = true;
                    beta_mean(condnan) = nan;
                    beta_var(condnan) = nan;
                    snr_mean(condnan) = nan;
                    snr_var(condnan) = nan;
                    velo_mean(condnan) = nan;
                    snr_instr_var(condnan) = nan;
                    
                    
                    hf = figure; hf.Units = 'centimeters'; hf.Position = [.5 2 25 10];
                    hf.Color = 'white'; hf.Visible = 'off';
                    
                    sp1 = subplot(321);
                    pcolor(data.time_3min,data.height/1000,real(log10(beta_mean))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 7.3 11 2.2]);
                    caxis([-7 -3]); colormap(sp1,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'Beta mean');
                    cb1 = colorbar; cb1.Label.String = 'm-1 sr-1'; ax1 = get(gca,'Position'); cb1.Units = 'centimeters';
                    cb1.Ticks = -7:-3; cb1.TickLabels = [repmat('10^{',length(cb1.Ticks(:)),1), ...
                        num2str(cb1.Ticks(:)) repmat('}',length(cb1.Ticks(:)),1)];
                    cb1.Position(3) = .25; cb1.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)');
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp2 = subplot(322);
                    pcolor(data.time_3min,data.height/1000,real(log10(beta_var))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 7.3 11 2.2]);
                    caxis([-16 -8]); colormap(sp2,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'Beta variance')
                    cb2 = colorbar; cb2.Label.String = '-'; ax1 = get(gca,'Position'); cb2.Units = 'centimeters';
                    cb2.Position(3) = .25; cb2.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb2.Ticks = -16:2:8; ylabel('Height (km)'); cb2.TickLabels = [repmat('10^{',length(cb2.Ticks(:)),1), ...
                        num2str(cb2.Ticks(:)) repmat('}',length(cb2.Ticks(:)),1)];
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp3 = subplot(323);
                    pcolor(data.time_3min,data.height/1000,10*real(log10(snr_mean-1))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 4.2 11 2.2]);
                    caxis([-30 10]); colormap(sp3,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'Signal mean');
                    cb3 = colorbar; cb3.Label.String = 'dB'; cb3.Ticks = -30:10:10; ax1 = get(gca,'Position'); cb3.Units = 'centimeters';
                    cb3.Position(3) = .25; cb3.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)')
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp4 = subplot(324);
                    pcolor(data.time_3min,data.height/1000,real(log10(snr_var))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 4.2 11 2.2]);
                    caxis([-7 1]); colormap(sp4,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'Signal variance')
                    cb4 = colorbar; cb4.Label.String = '-'; ax1 = get(gca,'Position'); cb4.Units = 'centimeters';
                    cb4.Position(3) = .25; cb4.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb4.Ticks = -7:2:1; cb4.TickLabels = [repmat('10^{',length(cb4.Ticks(:)),1), ...
                        num2str(cb4.Ticks(:)) repmat('}',length(cb4.Ticks(:)),1)];
                    ylabel('Height (km)')
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp5 = subplot(325);
                    pcolor(data.time_3min,data.height/1000,velo_mean'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[1 1.1 11 2.2]);
                    caxis([-2 2]); colormap(sp5,cmocean('balance')); text(0,3.35,'Velocity mean')
                    cb5 = colorbar; cb5.Label.String = 'm s-1'; ax1 = get(gca,'Position'); cb5.Units = 'centimeters';
                    cb5.Position(3) = .25; cb5.Position(1) = 10.3; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    cb5.Ticks = -2:.5:2;
                    ylabel('Height (km)'); xlabel('Time UTC')
                    set(gca,'Color',[.5 .5 .5])
                    
                    sp6 = subplot(326);
                    pcolor(data.time_3min,data.height/1000,real(log10(snr_instr_var))'); axis([0 24 0 3]); shading flat
                    set(gca,'Ytick',0:3,'XTick',0:3:24,'Units','centimeters','Position',[13.5 1.1 11 2.2]);
                    caxis([-7 1]); colormap(sp6,flipud(cbrewer('div','Spectral',256))); text(0,3.35,'Signal instrumental precision variance')
                    cb6 = colorbar; cb6.Label.String = '-'; ax1 = get(gca,'Position'); cb6.Units = 'centimeters';
                    cb6.Position(3) = .25;   cb6.Position(1) = 22.8; pause(.1); set(gca,'Position',ax1,'Units','centimeters');
                    ylabel('Height (km)'); xlabel('Time UTC'); cb6.Ticks = -7:2:1;
                    cb6.TickLabels = [repmat('10^{',length(cb6.Ticks(:)),1), ...
                        num2str(cb6.Ticks(:)) repmat('}',length(cb6.Ticks(:)),1)];
                    set(gca,'Color',[.5 .5 .5])
                    
                    if exist('typeof','var') == 1
                        [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode,typeof);
                        export_fig('-png',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                            '-%s-%s.png'], dir_out,num2str(DATE),site,measmode,typeof))
                    else
                        [dir_out,~] = getHALOfileList(site,DATE,processlev,measmode);
                        export_fig('-png','-m2',sprintf(['%s%s_%s_halo-doppler-lidar-' num2str(C.halo_unit_id) ...
                            '-%s.png'], dir_out,num2str(DATE),site,measmode))
                    end
                    close(hf)
                otherwise
                    continue
            end
    end
end
end

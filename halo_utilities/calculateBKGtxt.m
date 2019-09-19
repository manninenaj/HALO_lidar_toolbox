function [bkg_out, fit_out, bkg_times] = calculateBKGtxt(bkg_path,files_bkg,file_type,daten,n_range_gates,range_m)

irange=find(range_m<110,1,'first');
% find co and cross background files
dates=datestr(daten,'ddmmyy');
   
%%%filess=dir([bkg_path 'Background_' dates '*.txt']);


%switch file_type
 % case 'txt'
bkg_times = nan(length(files_bkg),1); % col1: time
if isempty(files_bkg)
    bkg = nan(1,n_range_gates);
else
    bkg = nan(length(files_bkg),n_range_gates);
end

switch file_type
    case 'txt'   
        bkg_times = nan(length(files_bkg),1); % col1: time
        for i = 1:length(files_bkg)
        %%%b_daten=datenum([filess(i).name(12:17) filess(i).name(19:24)],'ddmmyyHHMMSS');
            b_daten = datenum([files_bkg{i}(12:17) files_bkg{i}(19:24)],'ddmmyyHHMMSS');
            bkg_times(i,:) = b_daten;
    
            fn = [bkg_path files_bkg{i}];
            fid=fopen(fn,'r');
            bk=fscanf(fid,'%s');
            fclose(fid);
        
            dot_j=find(bk=='.');
            end_j=[1;transpose(dot_j+7)];
        
            bj=1;
            for jj=2:n_range_gates+1
                bkg(i,bj)=str2num(bk(end_j(jj-1):end_j(jj)-1));
                bj=bj+1;
            end
         end
  case 'nc' % blindly assume only one file, only ARM uses this..
    tmp = load_nc_struct([bkg_path,files_bkg{1}]);
    bkg_times = decimal2daten(tmp.time/3600,daten);
    bkg = tmp.background;
end

%% gapfilling

%if isempty(bkg) % no data for a day
%    bkg=nan(24,n_range_gates+1);
%    bkg(:,1)=(0:23)/24+daten;
%else
%    
%    bkg=[bkg_times bkg];
%    if str2num(datestr(bkg(1,1),'HH'))~=0
%        bkg=[[floor(bkg_times(1));bkg(:,1)] [bkg(1,2:end)*nan;bkg(:,2:end)]];
%    end
%    if str2num(datestr(bkg_times(end),'HH'))~=23
%        bkg=[[bkg(:,1); daten+23/24] [bkg(:,2:end);bkg(1,2:end)*nan]];
%    end
    
%    dd=(diff(bkg(:,1)-daten))*24;
%    for i=1:length(dd)
%        if dd(i)>1.5
%            new_t=trnspose((bkg(i,1)+1/24):1/24:(bkg(i,1)+floor(dd(i))/24));
%            bkg=[[bkg(:,1); new_t] [bkg(:,2:end);repmat(bkg(1,2:end)*nan,length(new_t),1)]];
%        end
%    end
    
%    if any(dd>1.5)
%        bkg=sortrows(bkg,1);
%    end
%end
   
bkg_raw=bkg;
%clear bkg;
%% to SNR
% bkg_snr=bkg_raw;
% bkg_snr(:,2:end)=nan;
fits_1=bkg_raw(:,1:3)*nan; % p(1) p(2) rmse
fits_2=bkg_raw(:,1:4)*nan; % p(1) p(2) p(3) rmse
fit_out = nan(length(bkg_raw(:,1)),n_range_gates);
bkg_out = nan(length(bkg_raw(:,1)),n_range_gates);
for i=1:length(bkg_raw(:,1));
    %b_temp=bkg_raw(i,2:end);
    %b_temp = b_temp(:);
    if ~isnan(bkg_raw(i,1)) % fit 1 and 2 order polynomial
        switch file_type
          case 'txt'
            x = irange:n_range_gates-1;
            y = bkg_raw(i,irange:n_range_gates-1);
            fitti_1 = my_robustfit(x(:),y(:));
            fitti_2 = my_robustfit([x(:) x(:).^2],y(:));
	    p1 = [fitti_1(2) fitti_1(1)];
            p2 = [fitti_2(3) fitti_2(2) fitti_2(1)];
            bkg_fitted_1 = polyval(p1, 1:n_range_gates);
            bkg_fitted_2 = polyval(p2, 1:n_range_gates);
            bkg_fitted_1 =  transpose(bkg_fitted_1(:));
            bkg_fitted_2 =  transpose(bkg_fitted_2(:)); 
            rmse_1 = sqrt(mean((bkg_raw(i,irange:n_range_gates)-bkg_fitted_1(irange:n_range_gates)).^2,2));
            rmse_2 = sqrt(mean((bkg_raw(i,irange:n_range_gates)-bkg_fitted_2(irange:n_range_gates)).^2,2));
          case 'nc'
            x = 1:length(bkg_raw(i,5:end)); x = transpose(x(:));
            fitti_1 = polyfit(x,bkg_raw(i,irange:end),1);
            fitti_2 = polyfit(x,bkg_raw(i,irange:end),2);
            bkg_fitted_1 = (1:n_range_gates)*fitti_1(1)+fitti_1(2);
            bkg_fitted_2 = ((1:n_range_gates).^2)*fitti_2(1)+(1:n_range_gates)*fitti_2(2)+fitti_2(3);
            bkg_fitted_1 =  bkg_fitted_1(:);
            bkg_fitted_2 =  bkg_fitted_2(:);
            rmse_1 = sqrt(mean((bkg_raw(i,irange:end)-bkg_fitted_1(irange:end)).^2,2));
            rmse_2 = sqrt(mean((bkg_raw(i,irange:end)-bkg_fitted_2(irange:end)).^2,2));
        end            
        if rmse_2<(0.9*rmse_1)
            fit_out(i,:) = bkg_fitted_2;
        else
            fit_out(i,:) = bkg_fitted_1;
        end
        
       switch file_type
         case 'txt'
           bkg_out(i,:) = bkg_raw(i,irange:n_range_gates);
         case 'nc'
           bkg_out(i,:) = bkg_raw(i,:);
       end
    end
end

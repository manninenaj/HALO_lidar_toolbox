function [bkg_out, fit_out, bkg_times] = calculateBKGtxt(bkg_path,daten,n_range_gates)

%% find co and cross background files

dates=datestr(daten,'ddmmyy');

filess=dir([bkg_path 'Background_' dates '*.txt']);

bkg_times=nan(length(filess),1); % col1: time
for i=1:length(filess)
    b_daten=datenum([filess(i).name(12:17) filess(i).name(19:24)],'ddmmyyHHMMSS');
    bkg_times(i,:)=b_daten;
end

%% read in backgrounds
bkg=nan(length(filess),n_range_gates);
for i=1:length(filess)
    fn=[bkg_path filess(i).name];
    fid=fopen(fn,'r');
    bk=fscanf(fid,'%s');
    fclose(fid);
    
    dot_i=find(bk=='.');
    end_i=[1;transpose(dot_i+7)];

    bi=1;
    %         for ii=2:length(end_i)
    for ii=2:n_range_gates+1
        bkg(i,bi)=str2num(bk(end_i(ii-1):end_i(ii)-1));
        bi=bi+1;
    end
end

%% gapfilling

if isempty(bkg) % no data for a day
    bkg=nan(24,n_range_gates+1);
    bkg(:,1)=(0:23)/24+daten;
else
    
    bkg=[bkg_times bkg];
    if str2num(datestr(bkg(1,1),'HH'))~=0
        bkg=[[floor(bkg_times(1));bkg(:,1)] [bkg(1,2:end)*nan;bkg(:,2:end)]];
    end
    if str2num(datestr(bkg_times(end),'HH'))~=23
        bkg=[[bkg(:,1); daten+23/24] [bkg(:,2:end);bkg(1,2:end)*nan]];
    end
    
    dd=(diff(bkg(:,1)-daten))*24;
    for i=1:length(dd)
        if dd(i)>1.5
            new_t=transpose((bkg(i,1)+1/24):1/24:(bkg(i,1)+floor(dd(i))/24));
            bkg=[[bkg(:,1); new_t] [bkg(:,2:end);repmat(bkg(1,2:end)*nan,length(new_t),1)]];
        end
    end
    
    if any(dd>1.5)
        bkg=sortrows(bkg,1);
    end
end
bkg_raw=bkg;
clear bkg;
%% to SNR
% bkg_snr=bkg_raw;
% bkg_snr(:,2:end)=nan;
fits_1=bkg_raw(:,1:3)*nan; % p(1) p(2) rmse
fits_2=bkg_raw(:,1:4)*nan; % p(1) p(2) p(3) rmse

fit_out = nan(length(bkg_raw(:,1)),n_range_gates);
bkg_out = nan(length(bkg_raw(:,1)),n_range_gates);
for i=1:length(bkg_raw(:,1));
    b_temp=bkg_raw(i,2:end);
    if ~isnan(b_temp(1)) % fit 1 and 2 order polynomial
        fitti_1 = polyfit((4:n_range_gates),b_temp(4:n_range_gates),1);
        bkg_fitted_1 = (1:n_range_gates)*fitti_1(1)+fitti_1(2);
        fitti_2=polyfit((4:n_range_gates),b_temp(4:n_range_gates),2);
        bkg_fitted_2=((1:n_range_gates).^2)*fitti_2(1)+(1:n_range_gates)*fitti_2(2)+fitti_2(3);
        rmse_1=sqrt(mean((b_temp(4:n_range_gates)-bkg_fitted_1(4:n_range_gates)).^2,2));
        rmse_2=sqrt(mean((b_temp(4:n_range_gates)-bkg_fitted_2(4:n_range_gates)).^2,2));
                    if rmse_2<(0.9*rmse_1)
                        fit_out(i,:) = bkg_fitted_2;
                    else
                        fit_out(i,:) = bkg_fitted_1;
                    end       
        bkg_out(i,:) = b_temp;
    end
end
end





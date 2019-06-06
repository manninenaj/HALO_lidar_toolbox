function droplet_bit = get_droplet_bit_peaks(height, beta)

% This function finds the bases of liquid clouds (droplet bit)

% Modified version of the original Cloudnet droplet bit function
% By Minttu Tuononen, Finnish Meteorological Institute
% 23.06.2017

% This fuction is using the findpeaks matlab function to find the maximum
% beta value and checks the base similarly as Cloudnet function

beta0 = beta;
beta0(find(~isfinite(beta0))) = 0;

droplet_bit = zeros(size(beta));
dheight = median(diff(height));
nrays = size(beta,1);
ngates = length(height);

threshold_beta = 2e-5;
min_jump = 10; % Factor of 20
jump_distance = 250; % m                               
jump_gates = ceil(jump_distance / dheight); %250
final_gate = ngates-jump_gates;
search_gates_below = ceil(100/dheight); %100 m
search_gates_above = ceil(300/dheight); %300 m
radar_search_gates_above = ceil(2000/dheight);
diff_factor = 0.25;
step_gates = 2;
min_temperature = 273.16 - 40; % No droplets colder than this

%return;

% fix to identify fog or strong precipitation:
%beta0(:,1)=0;
  
disp('Calculating location of liquid water cloud droplets')
for ii = 1:nrays

  start_gate = 2;
  profile = beta0(ii,:);
  
  % find peaks in the beta profile:
  % THINK ABOUT THIS: I might not want to add minpeakdistance
  % [mhighbeta,ihighbeta]=findpeaks(profile(:,:),'MinPeakProminence',2e-5,'MaxPeakWidth',15,'MinPeakDistance',30);
  
  % MinPeakProminence               - how large peak in magnitude
  % MinPeakDistance                 - how far away from another peak
  % WidthReference, 'halfheight'    - this measures the width based on the
  %                                   height of the peak (not prominence)
  % MaxPeakWidth                    - this sets the minimum width of the
  %                                   peak
  
  [mhighbeta,ihighbeta]=findpeaks(profile(:,:),'MinPeakProminence',2e-5,'MinPeakDistance',30,'WidthReference','halfheight','MaxPeakWidth',15);
 
      % Candidate liquid water
      % Find the base
      
      for i = 1:length(ihighbeta)
        %search_start = max(1, ihighbeta(i)-search_gates_below);
        if i == 1
          % this is original
          search_start = max(1, ihighbeta(i)-search_gates_below);
        else
          % in case of several peaks:
          % check that search starts either 250 m below or 
          % from previous peak depending on which is higher
          % two peaks cannot have same base
          search_start = max(1, max(ihighbeta(i)-search_gates_below, ihighbeta(i-1))); 
        end
        diff_profile = diff(profile(search_start:ihighbeta(i)));
        max_diff = max(diff_profile);
        
        % find if there is positive difference below
        % --> do not take into account as it belongs to a weak peak values
        
        % this is original cloudnet approach
        pot_base_liquid = min(find(diff_profile > max_diff*diff_factor))+search_start;
        
        
        % MIN IS A BIT PROBLEMATIC:
        % there may be additional, weaker peaks. 
        
        % check that the base is assosiated to the major peak:
        
        % ind neg tells if there are other peaks
        ind_neg = find(diff_profile<0);
        
      % jos base_liquid pienemmän piikin kohdalla --> korjaa ylöspäin
       if ~isempty(find(ind_neg == pot_base_liquid))
           ix = find(diff_profile > max_diff*diff_factor);
           base_liquid = min(ix(find(ix > pot_base_liquid)))+search_start;
       else 
           base_liquid = pot_base_liquid;
       end
        
        % check that beta value at base_liquid is smaller than peak value
       

       if profile(base_liquid)> mhighbeta(i)       
           if ~isempty(ind_neg)
            % start from last negative value --> add max(ind)
            base_liquid = min(find(diff_profile(max(ind_neg):end) > max_diff*diff_factor))+max(ind_neg)-1+search_start;
           end        
       end
   
  
  %         elseif profile(base_liquid)< mhighbeta(i) && ~isempty(ind_neg)
  %            base_liquid = min(find(diff_profile(max(ind_neg):end) > max_diff*diff_factor))+max(ind_neg)+search_start;
  %         end
  %       end
        
%         ind_neg = find(diff_profile<0);
%         if ~isempty(ind_neg)
%           % start from last negative value --> add max(ind)
%           base_liquid = min(find(diff_profile(max(ind_neg):end) > max_diff*diff_factor))+max(ind_neg)-2+search_start;
%         % else is original
%         else base_liquid = min(find(diff_profile > max_diff*diff_factor))+search_start;
%         end

          % Find the top
        % First see if profile goes to zero
%         top_liquid = min(find(~profile(ihighbeta(i)+[1:search_gates_above])))+ihighbeta(i)-1;
%           if isempty(top_liquid)
%             diff_profile = diff(profile(ihighbeta(i)+[0:search_gates_above]));
%             max_diff = max(-diff_profile);
%             top_liquid = max(find(-diff_profile > max_diff*diff_factor))+ihighbeta(i)-1;
%           end
      
        % Set to wet
        droplet_bit(ii, base_liquid) = 1;
      end
      
  %  end
  %  end
    % Previously we did this which seems very inefficient: repeated
    % analysis of the same pixels:
    %    start_gate+step_gates;
    % This seems more sensible:
    %start_gate = ihighbeta+1;
  %end
end

%droplet_bit(find(temperature < min_temperature)) = 0;
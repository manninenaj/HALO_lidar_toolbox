function [C,O] = getconfig(site,DATE)
%GETCONFIG reads the config file and outputs the site and date specific
%paremeters in struct format.
%
% Usage:
% C = getconfig(site,DATE)
% [C,fnames] = getconfig(site,DATE)
%
% Inputs:
% - site           string, name of the site, e.g. 'kuopio'
% - DATE           scalar, numerical date, e.g. 20171231
%
% Outputs:
% - C              struct, paramteres specific for site and data
% - O              cell array, 1st column contains names of different 
%                  measurement modes and the 2nd the type as given in the
%                  halo_config.txt file for the site and data
%
% Created 2017-09-31
% Antti Manninen
% Finnish Meteorological Institute
% antti.manninen(at)fmi.fi

% Check inputs
if nargin < 2
    error('''site'' and ''DATE'' are required inputs!')
end
if ~ischar(site)
    error('The first input ''site'' must be a string.')
end
if ~isnumeric(DATE) || length(num2str(DATE))~=8 
    error(['The second input ''DATE'' must be a numerical date in' ...
        ' YYYYMMDD format.'])
end

% Check that the config file can be found
if exist('halo_config.txt','file') == 2
    
    %%--- open, read, close ---%%
    fid = fopen('halo_config.txt');
    fspec = '%s';
    D = textscan(fid,fspec,'Delimiter','\n','HeaderLines',1);
    fclose(fid);
    D = D{:};
    
    % Sort parameter/value pairs into columns
    D1 = D(not(cellfun('isempty',D))); % rm empty lines
    i_keep = find(cellfun('isempty',strfind(D1,'SITE SPECIFIC')) & ~strncmp(D1,'%',1));
    D1 = D1(i_keep); % remove comment lines 
    split1 = regexp(D1,' = ','split'); % split after equals sign
    param_val_pairs = cell(length(split1),2);
    for i=1:length(split1)
        param_val_pairs{i,1} = split1{i,1}{1,1};
        param_val_pairs{i,2} = split1{i,1}{1,2};
    end
    
    % Check that the site can be found from the halo_config.txt file
    if any(strcmp(param_val_pairs(:,2),site))
    
        %%--- Get default parameters ---%%
        % Set limits
        iend_def = find(strcmp(param_val_pairs(:,1),'site'),1,'first')-1;
        
        % Get names and values and allocate them into a struct variable
        def_param_names = param_val_pairs(1:iend_def,1);
        def_param_values = param_val_pairs(1:iend_def,2);
        for i = 1:length(def_param_names)
            if ~isnan(str2double(def_param_values{i}))
                % Convert numeric parameters into double precision
                C.(def_param_names{i}) = ...
                    str2double(def_param_values{i});
%                 C.default.(def_param_names{i}) = ...
%                     str2double(def_param_values{i});
            else
                % Text as is
                C.(def_param_names{i}) = ...
                    def_param_values{i};
%                 C.default.(def_param_names{i}) = ...
%                     def_param_values{i};
            end
        end
        %%--- Get site specific parameters ---%%
        % Set limits
        ibegin = find(strcmp(param_val_pairs(:,2),site),1,'first');
        iend = ibegin + find(strcmp(param_val_pairs(ibegin+1:end,1),'site'),1,'first')-1;

        if isempty(iend), iend = length(param_val_pairs); end % in case last
    
        % Get names and values
        spcfc_param_names = param_val_pairs(ibegin:iend,1);
        spcfc_param_values = param_val_pairs(ibegin:iend,2);
        
        % Check parameters_valid_from_including for the time period
        iperiods = find(strcmp(spcfc_param_names,...
            'parameters_valid_from_including'));

        if any(datenum(num2str(DATE),'yyyymmdd') >= ...
            datenum(spcfc_param_values(iperiods),'yyyymmdd'))

        % Read parameter-value pairs for the appropriate date
        valid_from_date = '0';
        for j = 1:length(spcfc_param_values)
          % Update valid time period from 'parameters_valid_from_including'
          if strcmp(spcfc_param_names{j},'parameters_valid_from_including')
            valid_from_date = spcfc_param_values{j};
          end

          if datenum(num2str(DATE),'yyyymmdd') >= datenum(valid_from_date,'yyyymmdd')
            % Check if value is a number or a string
            if ~isnan(str2double(spcfc_param_values{j}))
              C.(spcfc_param_names{j}) = str2double(spcfc_param_values{j});
            else
              C.(spcfc_param_names{j}) = spcfc_param_values{j};
            end
          end
        end       





%	  for i = 1:length(iperiods)
%              %
%
%              if datenum(num2str(DATE),'yyyymmdd') >= ...
%	          datenum(spcfc_param_values(iperiods(i)),'yyyymmdd')
%		  for j = 1:length(spcfc_param_values)
%                      % number
%                      if ~isnan(str2double(spcfc_param_values{j}))
%	                  C.(spcfc_param_names{j}) = ...
%		              str2double(spcfc_param_values{j});
%                      % string
%                      else
%                          C.(spcfc_param_names{j}) = ...
%                              spcfc_param_values{j};
%                      end
%		  end
%              end
%	  end	    		    
        else
            error(['DATE = %d < the earliest valid date specified in' ...
                  ' the halo_config.txt\nfor'' %s'' site.' ...
                  ' Please check ''parameters_valid_from_including''' ...
                  '\nfrom halo_config.txt.'],DATE,site)
        end
    else
        error(['''%s'' is not specified in halo_config.txt.\nCheck' ...
            ' spelling or add your site to the halo_config.txt.'],site)
    end

else
    error('Please check that halo_config.txt is in your path.')
end

% Output the names of different measurements
a0 = fieldnames(C);
% a0 = fieldnames(C.site_specific);
a1 = a0(~cellfun('isempty',strfind(a0,'dir_original')));
a2 = cellfun(@(x) x(14:end),a1,'UniformOutput',false);
b1 = cellfun(@(x,y) x(1:y-1),a2,strfind(a2,'_'),'UniformOutput',false);
c1 = cellfun(@(x,y) x(y+1:end),a2,strfind(a2,'_'),'UniformOutput',false);
O(:,1) = b1; O(:,2) = c1;


% Order fields into alphabetical order
C = orderfields(C);
% C.default = orderfields(C.default);
% C.site_specific= orderfields(C.site_specific);

end

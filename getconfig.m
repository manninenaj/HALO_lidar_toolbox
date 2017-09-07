function C = getconfig(site,DATE)
%GETCONFIG reads the config file and outputs the site and date specific
%paremeters in struct format.

% Check inputs, TBD: improve
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
        def_param_values = str2double(param_val_pairs(1:iend_def,2));
        for i = 1:length(def_param_names)
            C.(def_param_names{i}) = def_param_values(i);
        end
        
        %%--- Get site specific parameters ---%%
        % Set limits
        ibegin = find(strcmp(param_val_pairs(:,2),site))+1;
        iend = find(strcmp(param_val_pairs(ibegin+1:end,1),'site')) + ...
            ibegin-1;
        if isempty(iend), iend = length(param_val_pairs); end % if last
        
        % Get names and values
        spcfc_param_names = param_val_pairs(ibegin:iend,1);
        spcfc_param_values = param_val_pairs(ibegin:iend,2);
        
        % Check parameters_valid_from_including for the time period
        iperiods = find(strcmp(spcfc_param_names,...
            'parameters_valid_from_including'));
        
        if datenum(num2str(DATE),'yyyymmdd') >= ...
                datenum(spcfc_param_values(iperiods(1)),'yyyymmdd')
            
            iend2 = iperiods(find(datenum(num2str(DATE),'yyyymmdd') < ...
                datenum(spcfc_param_values(iperiods),'yyyymmdd')...
                ,1,'last'))-1;
            
            % if last period
            if isempty(iend2), iend2 = length(spcfc_param_values); end
            
            % Allocate parameter/value pairs into a struct variable,
            % overwrite default parameters if they are different for  
            % the site and the time period
            for i = 1:iend2
                if ~isnan(str2double(spcfc_param_values{i}))
                    % Convert numeric parameters into double precision
                    C.(spcfc_param_names{i}) = ...
                        str2double(spcfc_param_values{i});
                else
                    % Text as is
                    C.(spcfc_param_names{i}) = ...
                        spcfc_param_values{i};
                end
            end
        else
            error(['DATE = %d < the earliest valid date specified in' ...
                ' the halo_config.txt\nfor'' %s'' site.' ...
                ' Please check ''parameters_valid_from_including''' ...
                '\nfrom halo_config.txt.'],DATE,site)
        end
    else
        error(['Site ''%s'' not recognized as a site specified in' ...
            ' halo_config.txt.\nPlease check spelling of ''site'''...
            ' or specify your site in halo_config.txt.'],site)
    end
else
    error('Please check that halo_config.txt is in your path.')
end
end


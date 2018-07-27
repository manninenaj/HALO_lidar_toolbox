function []=write_nc_silent(nc_file, dimensions, data, attributes)

disp(['Writing ' nc_file]);
f = netcdf.create(nc_file, 'clobber');
% mode = netcdf.getConstant('NETCDF4');
% mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
% f = netcdf.create(nc_file, mode);

% Define dimensions
dimension_names = fieldnames(dimensions);
for ii = 1:length(dimension_names)
%     disp(['   Adding dimension ' dimension_names{ii}]);
    netcdf.defDim(f,dimension_names{ii},getfield(dimensions, dimension_names{ii}));
    %   f(dimension_names{ii}) = getfield(dimensions, dimension_names{ii});
end
% scalar_dim=netcdf.defDim(f,'scalar_dimension',1);

% Set global attributes
varid = netcdf.getConstant('GLOBAL');
attribute_names = fieldnames(attributes.global);
for ii = 1:length(attribute_names)
    %     disp(['   Adding global attribute ' attribute_names{ii}]);
    value = getfield(attributes.global, attribute_names{ii});
    [nc_class, straightjacket] = class2ncclass(value);
    if ~isempty(nc_class)
        if strcmp(nc_class,'float')
            netcdf.putAtt(f,varid,attribute_names{ii},single(value))
        else
            netcdf.putAtt(f,varid,attribute_names{ii},value)
        end
        %     eval(['f.' attribute_names{ii} ' = nc' nc_class '(' straightjacket '(value));']);
    else
        warning(['Global attribute ' attribute_names{ii} ' omitted due to incompatible type']);
    end
end

% Set variable information
variable_names = fieldnames(data);
for ii = 1:length(variable_names)
    disp(['   Adding variable ' variable_names{ii}]);
    value = getfield(data, variable_names{ii});
    varattributes = getfield(attributes, variable_names{ii});
    vardimensions = getfield(varattributes, 'dimensions');
    nc_class = class2ncclass(value);
    if ~isempty(nc_class)
        if isfield(varattributes,'missing_value')
            missing_value = varattributes.missing_value;
        end
        %     varinfo = {nc_class};
%         dimids=scalar_dim;
        dimids=[];
        for jj = 1:length(vardimensions)
            dimids(jj) = netcdf.inqDimID(f, vardimensions{jj});
            %       varinfo{jj+1} = vardimensions{jj};
        end
        try
            %         varid=netcdf.defVar(f,variable_names{ii},nc_class,dimids);
            varid=netcdf.defVar(f,variable_names{ii},nc_class,dimids(end:-1:1));
        catch
            % if length(dimids)==2
            keyboard
            % end
        end
        %     f{variable_names{ii}} = varinfo;
        
        % Set attributes
        attribute_names = fieldnames(varattributes);
        for jj = 1:length(attribute_names)
            if strcmp(attribute_names{jj}, 'dimensions')
                continue;
            end
%             disp(['      adding attribute ' variable_names{ii} '.' attribute_names{jj}]);
            value = getfield(varattributes, attribute_names{jj});
            
%             netcdf.putAtt(f,varid,attribute_names{jj},value);
            
            [attr_nc_class, straightjacket] = class2ncclass(value);
            if ~isempty(attr_nc_class)
                if strcmp(attribute_names{jj},'FillValue_')
                    % netcdf.defVarFill(f,varid,false,single(value));
                    netcdf.putAtt(f,varid,'_FillValue',single(value));
                else
                    if strcmp(attribute_names{jj},'missing_value') | strcmp(attribute_names{jj},'plot_range')
                        attr_nc_class = nc_class;
                    end
                    if strcmp(attr_nc_class,'float')
                        netcdf.putAtt(f,varid,attribute_names{jj},single(value));
                    else
                        netcdf.putAtt(f,varid,attribute_names{jj},value);
                    end
                end
%                 eval(['f{variable_names{ii}}.' attribute_names{jj} ' = nc' attr_nc_class '(' straightjacket '(value));']);
            else
                warning(['Attribute ' variable_names{ii} '.' attribute_names{jj} ' omitted due to incompatible type']);
            end
        end
    
        % Fill variables
        netcdf.endDef(f)
        data_value = getfield(data, variable_names{ii});
        data_varattributes = getfield(attributes, variable_names{ii});
        nc_class = class2ncclass(data_value);
%         disp(['   Filling variable ' variable_names{ii} ' (' nc_class ')']);
        if ~isempty(nc_class)
            if isfield(data_varattributes,'missing_value') & strcmp(nc_class,'float')
                missing_value = data_varattributes.missing_value;
                data_value(~isfinite(data_value)) = missing_value;
            end
%             try
                if isempty(dimids)
                    if strcmp(nc_class,'float')
                        netcdf.putVar(f,varid,single(data_value'))
                    else
                        netcdf.putVar(f,varid,data_value')
                    end
                else
                    data_start=0;
                    data_len=0;
                    for ij=1:length(dimids)
                        data_start(ij)=0;
                        data_len(ij)=size(data_value,ij);
                    end
                    if strcmp(nc_class,'float')
                        if length(dimids)==3
                            netcdf.putVar(f,varid,data_start(end:-1:1),data_len(end:-1:1),single(permute(data_value,[3 2 1])))
                        else
                            netcdf.putVar(f,varid,data_start(end:-1:1),data_len(end:-1:1),single(data_value'))
                        end
                    else
                        if length(dimids)==3
                            netcdf.putVar(f,varid,data_start(end:-1:1),data_len(end:-1:1),permute(data_value,[3 2 1]))
                        else
                            netcdf.putVar(f,varid,data_start(end:-1:1),data_len(end:-1:1),data_value')
                        end
                    end
                end
%             catch e
%                 keyboard
%             end
%             f{variable_names{ii}}(:) = value;
        else
            warning(['Variable ' variable_names{ii} ' ommitted due to incompatible type']);
        end
        netcdf.reDef(f)
        
        
    else
        warning(['Variable ' variable_names{ii} ' ommitted due to incompatible type']);
    end
end

% % Fill variables
% variable_names = fieldnames(data);
% for ii = 1:length(variable_names)
%     value = getfield(data, variable_names{ii});
%     varattributes = getfield(attributes, variable_names{ii});
%     nc_class = class2ncclass(value);
%     disp(['   Filling variable ' variable_names{ii} ' (' nc_class ')']);
%     if ~isempty(nc_class)
%         if isfield(varattributes,'missing_value') & strcmp(nc_class,'float')
%             missing_value = varattributes.missing_value;
%             value(find(isnan(value))) = missing_value;
%         end
%         f{variable_names{ii}}(:) = value;
%     else
%         warning(['Variable ' variable_names{ii} ' ommitted due to incompatible type']);
%     end
% end

netcdf.close(f)

function [nc_class, straightjacket]  = class2ncclass(variable)
matclass = class(variable);
straightjacket = '';
if strcmp(matclass, 'double')
    nc_class = 'float'; % Only ever output floats...
elseif strcmp(matclass, 'single')
    nc_class = 'float';
elseif strcmp(matclass, 'char')
    nc_class = 'char';
elseif strcmp(matclass, 'int8')
    nc_class = 'byte';
    straightjacket = 'double';
elseif strcmp(matclass, 'int16')
    nc_class = 'short';
    straightjacket = 'double';
elseif strcmp(matclass, 'int32')
    nc_class = 'int';
else
    warning(['Unable to save matlab variables of type ' matclass ' in NetCDF']);
    nc_class = '';
end


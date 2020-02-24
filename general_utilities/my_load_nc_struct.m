function [data,att] = my_load_nc_struct(file_name, var_names)

data = [];
att = [];
if nargout < 1; error('No outputs..there must be outputs!'); end
info = ncinfo(file_name);
if nargin < 2
var_names = {info.Variables.Name};
end

for i = 1:length(var_names)
    data.(var_names{i}) = ncread(file_name, var_names{i});

    % Tranpose into desired order (time, height, ...)
    len_d = length(info.Variables(i).Dimensions); 
    if len_d == 2
        dim_names = {info.Variables(i).Dimensions.Name};
        itime = strmatch('time',dim_names);
        irange = strmatch('range',dim_names);
        if isempty(irange)
	    irange = strmatch('height',dim_names);
        end
        if isempty(irange)
	    irange = strmatch('height_agl',dim_names);
        end
        if isempty(irange)
            error('Unknown dimension name for 2nd dimension (expecting range, height, height_agl)!')
        end

        data.(var_names{i}) = permute(data.(var_names{i}), [itime, irange]);
    elseif len_d > 2
        error('Only works with 1-D or 2-D arrays!')
    end

    att_names = {info.Variables(i).Attributes.Name};
    att_vals = {info.Variables(i).Attributes.Value};
    for j = 1:length(att_names)
        att.(var_names{i}).(att_names{j}) = att_vals{j};
    end
end

end


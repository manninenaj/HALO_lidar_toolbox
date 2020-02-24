function [data,att] = my_load_nc_struct(file_name, var_names)

data = [];
att = [];
if nargout < 1; error('No outputs!'); end
info = ncinfo(file_name);
if nargin < 2
    var_names = {info.Variables.Name}';
end

for i = 1:length(var_names)
    data.(var_names{i}) = ncread(file_name, var_names{i});
    att_names = {info.Variables(i).Attributes.Name}';
    att_vals = {info.Variables(i).Attributes.Value}';
    for j = 1:length({info.Variables(i).Attributes.Name}')
        att.(var_names{i}).(att_names{j}) = att_vals{j};
    end
end

end


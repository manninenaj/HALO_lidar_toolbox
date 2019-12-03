function attribute = create_attributes(dimensions, long_name, units, missing_value, comment, plot)
attribute.dimensions = dimensions;
attribute.long_name = long_name;
if iscell(units)
    attribute.units = units{1};
    attribute.units_html = units{2};
else
    attribute.units = units;
end

if nargin > 3
    if ~isempty(missing_value)
        attribute.missing_value = missing_value;
        attribute.FillValue_ = missing_value;
    end
    if nargin > 4
        attribute.comment = comment;
        if nargin > 5
            if iscell(plot)
                attribute.plot_range = plot{1};
                attribute.plot_scale = plot{2};
            else
                attribute.plot_range = plot;
            end
        end
    end
end
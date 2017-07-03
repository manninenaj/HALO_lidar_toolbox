function [bad,wasnan,varargout]=statremovenan(varargin)
%STATREMOVENAN Remove NaN values from inputs



% expect three inputs (for now)!
y = varargin{1};
x = varargin{2};
priorw = varargin{3};

y_out = y;
x_out = x;
priorw_out = priorw;

bad = 0;
wasnan = isnan(y);
index = find(wasnan==1);
y_out(index) = [];
x_out(index) = [];
priorw_out(index) = [];
varargout{1} = y_out;
varargout{2} = x_out;
varargout{3} = priorw_out;

if size(y_out) ~= size(x_out)
  bad = 2;
end

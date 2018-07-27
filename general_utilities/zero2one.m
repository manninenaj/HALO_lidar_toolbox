function [ x_norm ] = zero2one(x)
%ZERO2ONE normalizes the values of the given array between 0 and 1
%
% 2015-11-24
% Antti Manninen
% antti.j.manninen@helsinki.fi
x_norm = (x - nanmin(x(:))) ./ (nanmax(x(:)) - nanmin(x(:)));
end


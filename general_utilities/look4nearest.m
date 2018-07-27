function [ab,ib] = look4nearest(a,b)
%look4nearest function looks for nearest values of vector 'a' in vector 'b'
%
% Inputs:
% - a       numeric vector, "look for"
% - b       numeric vector, "look from"
%
% Outputs:
% - ab      numeric vector, replacements for values in 'a' by nearest values in 'b'
% - ib      numeric vector, indices of 'b' nearest to the values in 'a' 
%
% Usage:
% ab = look4nearest(a,b)
% [~,ib] = look4nearest(a,b)
% [ab,ib] = look4nearest(a,b)
%
% Created 2017-12-14
% Antti Manninen
% University of Helsinki, Finland
% antti.j.manninen(at)helsinki.fi

% Check number of inputs
if nargin == 2
    % Check are the inputs vectors
    if isnumeric(a) && isnumeric(b) && isvector(a) && isvector(b)
        
        % Convert to row vectors
        a = a(:); a = reshape(a,1,length(a));
        b = b(:); b = reshape(b,1,length(b));
        
        [c,p] = sort(b);
        [~,ic] = histc(a,[-inf,(c(1:end-1)+c(2:end))/2,inf]);
        ib = p(ic);
        ab = b(ib);
        
    else
        disp('Inputs should be vectors!')
    end
else
    disp('Number of inputs should be 2 in look4nearest function!')
end
ab = ab(:);
ib = ib(:);
end


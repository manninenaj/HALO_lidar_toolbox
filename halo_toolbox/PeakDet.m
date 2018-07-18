function [MaxTab, MinTab]= PeakDet(Vect, Delta, Ex)
% PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% sendtome@billauer.co.il
% This function is released to the public domain; Any use is allowed.

MaxTab = [];
MinTab = [];

Vect = Vect(:); % Just in case this wasn't a proper vector

% Check inputs
if nargin < 3
    Ex = (1:length(Vect))';
else
    Ex = Ex(:);
    if length(Vect)~= length(Ex)
        error('Input vectors v and x must have same length');
    end
end

if (length(Delta(:)))>1
    error('Input argument DELTA must be a scalar');
end

if Delta <= 0
    error('Input argument DELTA must be positive');
end

% Prepare data
mn = Inf;
mx = -Inf;
mnpos = nan;
mxpos = nan;

lookformax = 1;

for i=1:length(Vect)
    this = Vect(i);
    if this > mx
        mx = this;
        mxpos = Ex(i);
    end
    if this < mn
        mn = this;
        mnpos = Ex(i);
    end
    
    if lookformax
        if this < mx-Delta
            MaxTab = [MaxTab ; mxpos mx];
            mn = this; mnpos = Ex(i);
            lookformax = 0;
        end
    else
        if this > mn+Delta
            MinTab = [MinTab ; mnpos mn];
            mx = this; mxpos = Ex(i);
            lookformax = 1;
        end
    end
end
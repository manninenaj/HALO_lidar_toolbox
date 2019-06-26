function [approx_coeff, detail_coeff] = stationaryWaveletDecomposition(x, w_level, w_type)
%stationaryWaveletTransform performs a stationary multilevel wavelet 
% decomposition by convolving the input signal with a selected filter
% function.
%
%Usage:
% [approx_coeff, detail_coeff] = my_swt(x,w_level,wavelet_type);
%
%Inputs:
% - x               vector, input signal to be denoised
% - w_level         scalar, wavelet decomposition level 
% - w_type          string, wavelet type ('haar', 'sym8')
%
%Outputs:
% - approx_coeff    vector, approximation coefficients for each level
% - detail_coeff    vector, detail coefficients for each level
%
% Created 2019-02-06
% Antti Manninen
% antti.j.manninen(at)helsinki.fi
% INAR / Physics
% University of Helsinki, Finland

% Use row vector.
x = transpose(x(:));
x_len = length(x);
% Check that the length of x is divisible by 2^(wavelet
% decomposition level)
if rem(x_len,2^w_level)>0
    warning('Input x should have the length, which is divisible by 2^(w_level). Check zero padding.');
    return
end

% Get decomposition filters.
switch w_type
    case {'haar','sym1'}
        low_pass = [ 1./sqrt(2) 1./sqrt(2)];
        high_pass = [-1./sqrt(2) 1./sqrt(2)];
    % Symlets
    case 'sym2'
        low_pass = [-0.1294, 0.2241, 0.8365, 0.4830];
        high_pass = [-0.4830, 0.8365, -0.2241, -0.1294];
    case 'sym3'
        low_pass = [0.0352, -0.0854, -0.1350, 0.4599, 0.8069, 0.3327];
        high_pass = [-0.3327, 0.8069, -0.4599, -0.1350, 0.0854, 0.0352];
    case 'sym4'
        low_pass = [-0.0758, -0.0296, 0.4976, 0.8037, 0.2979, -0.0992, -0.0126, 0.0322];
        high_pass = [-0.0322, -0.0126, 0.0992, 0.2979, -0.8037, 0.4976, 0.0296, -0.0758];
    case 'sym5'
        low_pass = [0.0273, 0.0295, -0.0391, 0.1994, 0.7234, 0.6340, 0.0166, -0.1753, -0.0211, 0.0195];
        high_pass = [-0.0195, -0.0211, 0.1753, 0.0166, -0.6340, 0.7234, -0.1994, -0.0391, -0.0295, 0.0273];
    case 'sym6'
        low_pass = [0.0154, 0.0035, -0.1180, -0.0483, 0.4911, 0.7876, 0.3379, -0.0726, -0.0211, 0.0447, 0.0018, -0.0078];
        high_pass = [0.0078, 0.0018, -0.0447, -0.0211, 0.0726, 0.3379, -0.7876, 0.4911, 0.0483, -0.1180, -0.0035, 0.0154];
    case 'sym7'
        low_pass = [0.0027, -0.0010, -0.0126, 0.0305, 0.0679, -0.0496, 0.0174, 0.5361, 0.7678, 0.2886, -0.1400, -0.1078, 0.0040, 0.0103];
        high_pass = [-0.0103, 0.0040, 0.1078, -0.1400, -0.2886, 0.7678, -0.5361, 0.0174, 0.0496, 0.0679, -0.0305, -0.0126, 0.0010, 0.0027];
    case 'sym8'
        low_pass = [-0.0034, -0.0005, 0.0317, 0.0076, -0.1433, -0.0613, 0.4814, 0.7772, 0.3644, -0.0519, -0.0272, 0.0491, 0.0038, -0.0150, -0.0003, 0.0019];
        high_pass = [-0.0019, -0.0003, 0.0150, 0.0038, -0.0491, -0.0272, 0.0519, 0.3644, -0.7772, 0.4814, 0.0613, -0.1433, -0.0076, 0.0317, 0.0005, -0.0034];
    % Daubechies
    case 'db4'
        low_pass = [-0.0106, 0.0329, 0.0308, -0.1870, -0.0280, 0.6309, 0.7148, 0.2304];
        high_pass = [-0.2304, 0.7148, -0.6309, -0.0280, 0.1870, 0.0308, -0.0329, -0.0106];
    case 'db5'
        low_pass = [0.0033, -0.0126, -0.0062, 0.0776, -0.0322, -0.2423, 0.1384, 0.7243, 0.6038, 0.1601];
        high_pass = [-0.1601, 0.6038, -0.7243, 0.1384, 0.2423, -0.0322, -0.0776, -0.0062, 0.0126, 0.0033];
    case 'db6'
        low_pass = [-0.0011, 0.0048, 0.0006, -0.0316, 0.0275, 0.0975, -0.1298, -0.2263, 0.3153, 0.7511, 0.4946, 0.1115];
        high_pass = [-0.1115, 0.4946, -0.7511, 0.3153, 0.2263, -0.1298, -0.0975, 0.0275, 0.0316, 0.0006, -0.0048, -0.0011];
    case 'db7'
        low_pass = [0.0004, -0.0018, 0.0004, 0.0126, -0.0166, -0.0380, 0.0806, 0.0713, -0.2240, -0.1439, 0.4698, 0.7291, 0.3965, 0.0779];
        high_pass = [-0.0779, 0.3965, -0.7291, 0.4698, 0.1439, -0.2240, -0.0713, 0.0806, 0.0380, -0.0166, -0.0126, 0.0004, 0.0018, 0.0004];
    case 'db8'
        low_pass = [-0.0001, 0.0007, -0.0004, -0.0049, 0.0087, 0.0140, -0.0441, -0.0174, 0.1287, 0.0005, -0.2840, -0.0158, 0.5854, 0.6756, 0.3129, 0.0544];
        high_pass = [-0.0544, 0.3129, -0.6756, 0.5854, 0.0158, -0.2840, -0.0005, 0.1287, 0.0174, -0.0441, -0.0140, 0.0087, 0.0049, -0.0004, -0.0007, -0.0001];
    otherwise
        warning('Only ''haar'', ''sym[1...8]'', or db[4...8] wavelets are supported.');
        return
end

% Compute stationary wavelet coefficients.
approx_coeff = zeros(w_level,x_len);
detail_coeff = zeros(w_level,x_len);

for k = 1:w_level
    
    % Extension
    lf = length(low_pass); % length of filter
    
    % Discrete Wavelet Transform mode is periodisation
    x = extendPeriodDWT(x,lf/2);
    
    % Decomposition
    detail_coeff(k,:) = extractVector(conv2(x(:)',high_pass(:)',...
        'full'),x_len,lf+1);
    approx_coeff(k,:) = extractVector(conv2(x(:)',low_pass(:)',...
        'full'),x_len,lf+1);
    
    % Dyadic upsampling of filters
    tmp = zeros(1,2.*lf);
    tmp(1:2:2 * lf) = low_pass;
    low_pass = tmp;
    tmp = zeros(1,2.*lf);
    tmp(1:2:2 * lf) = high_pass;
    high_pass = tmp;
    
    % Update x
    x = approx_coeff(k,:);
end
end

%-- Subfunction - extractVector
function y = extractVector(x, len, start)
%extractVector extracts a vector from within a larger vector

y = x;
finish = start + len - 1;
y = y(start:finish);

end

%-- Subfunction - extendPeriodDWT
function x = extendPeriodDWT(x,lf)
%extendPeriodDWT extends the DWT using periodisation

length_of_x = length(x);
if rem(length_of_x,2)
    x(length_of_x+1) = x(length_of_x);
    length_of_x = length_of_x+1;
end

I = [length_of_x-lf+1:length_of_x , 1:length_of_x , 1:lf];
if length_of_x<lf
    I = mod(I,length_of_x);
    I(I==0) = length_of_x;
end

x = x(I);

end


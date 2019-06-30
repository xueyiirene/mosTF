function [ Iout] = NormalizeImage(FIXED )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
% Normalize FIXED image

% Get linear indices to finite valued data
finiteIdx = isfinite(FIXED(:));
% 
% Replace NaN values with 0
FIXED(isnan(FIXED)) = 0;

% Replace Inf values with 1
FIXED(FIXED==Inf) = 1;

% Replace -Inf values with 0
FIXED(FIXED==-Inf) = 0;

% Normalize input data to range in [0,1].
FIXEDmin = min(FIXED(:));
FIXEDmax = max(FIXED(:));
if isequal(FIXEDmax,FIXEDmin)
    FIXED = 0*FIXED;
else
    FIXED(finiteIdx) = (FIXED(finiteIdx) - FIXEDmin) ./ (FIXEDmax - FIXEDmin);
end
Iout=FIXED;
end


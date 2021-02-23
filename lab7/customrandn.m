function randValueN = customrandn(m,s)
% Generates a normally distributed pseudo-random number
% R = customrandn(M,S)
% Generates a normally distributed pseudo-random number R.
% M is the mean of the normal distribution.
% S is the standard deviation of the normal distribution.
% Custom normal PDF N(x;mean,std) of Y = std*X+mean, where X has a PDF N(x;0,1).

% Alejandro Reyes, Feb 2021

randValueN = s*randn() + m;
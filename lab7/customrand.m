function randValue = customrand(a,b)
% Generates a uniformly distributed pseudo-random number
% R = customrand(A,B)
% Generates a pseudo-random number R between the open interval(A,B).
% A, B are the boundaries of the uniform probability distribution function.
% Custom uniform PDF U(x;a,b) of Y = (b-a)X+a, where X has a PDF U(x;0,1)

% Alejandro Reyes, Feb 2021

randValue = (b-a)*rand() + a;
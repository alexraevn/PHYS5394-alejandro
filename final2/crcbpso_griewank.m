function [fitVal,varargout] = crcbpso_griewank(xVec,params)
%A Griewank benchmark test function for CRCBPSO
%F = CRCBPSOTESTFUNC(X,P)
%Compute the Rastrigin fitness function for each row of X.  The fitness
%values are returned in F. X is standardized, that is 0<=X(i,j)<=1. P has
%two arrays P.rmin and P.rmax that are used to convert X(i,j) internally to
%actual coordinate values before computing fitness: X(:,j) ->
%X(:,j)*(rmax(j)-rmin(j))+rmin(j). 
%
%For standardized coordinates, F = infty if a point X(i,:) falls
%outside the hypercube defined by 0<=X(i,j)<=1.
%
%[F,R] =  CRCBPSOTESTFUNC(X,P)
%returns the real coordinates in R. 
%
%[F,R,Xp] = CRCBPSOTESTFUNC(X,P)
%Returns the standardized coordinates in Xp. This option is to be used when
%there are special boundary conditions (such as wrapping of angular
%coordinates) that are better handled by the fitness function itself.

% Based on crcbpsotestfunc.m by Soumya Mohanty
% Alejandro Reyes
% April 2021
%==========================================================================

%rows: points
%columns: coordinates of a point
[nrows,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nrows,1);
validPts = ones(nrows,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
%Convert valid points to actual locations
xVec(validPts,:) = s2rv(xVec(validPts,:),params);


for lpc = 1:nrows
    if validPts(lpc)
        % define coordinate row x
        x = xVec(lpc,:);
        % index vector to use inside the cosine of Griewank function
        ind = 1:length(x(1,:));
        % define the Generalized Griewank benchmark fitness function
        fitVal(lpc) = sum(x.^2)/4000-prod(cos(x./ind))+1;
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
    if nargout > 2
        varargout{2} = r2sv(xVec,params);
    end
end







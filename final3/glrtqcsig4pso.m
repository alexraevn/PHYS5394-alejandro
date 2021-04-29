function [fitVal,varargout] = glrtqcsig4pso(xVec,params)
%Fitness function for quadratic chirp regression
%F = CRCBQCFITFUNC(X,P)
%Compute the fitness function (log-likelihood ratio for colored noise maximized
%over the ampliture parameter) for data containing the quadratic chirp signal, X.
%The fitness values are returned in F. X is standardized, that is 0<=X(i,j)<=1.
%The fields P.rmin and P.rmax  are used to convert X(i,j) internally before
%computing the fitness:
%X(:,j) -> X(:,j)*(rmax(j)-rmin(j))+rmin(j).
%The fields P.dataY and P.dataX are used to transport the data and its
%time stamps. The fields P.dataXSq and P.dataXCb contain the timestamps
%squared and cubed respectively.
%The field P.psdVec is the PSD vector for the colored noise
%
%[F,R] = CRCBQCFITFUNC(X,P)
%returns the quadratic chirp coefficients corresponding to the rows of X in R. 
%
%[F,R,S] = CRCBQCFITFUNC(X,P)
%Returns the quadratic chirp signals corresponding to the rows of X in S.

%Alejandro Reyes
%Adapted from SDM/SDMBIGDAT19/CRCBQCFITFUNC
%==========================================================================

%rows: points
%columns: coordinates of a point
[nVecs,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nVecs,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
xVec(validPts,:) = s2rv(xVec(validPts,:),params);

for lpc = 1:nVecs
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        fitVal(lpc) = ssrqc(x, params);
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
end

%LLR after maximizing over amplitude parameter
function ssrVal = ssrqc(x,params)
%Generate normalized quadratic chirp
phaseVec = x(1)*params.dataX + x(2)*params.dataXSq + x(3)*params.dataXCb;
qc = sin(2*pi*phaseVec);
% Inner product of signal with itself
normSigSqrd = innerprodpsd(qc,qc,params.samplFreq,params.psdVec);
% Normalization factor
normFac = params.snr/sqrt(normSigSqrd);
% Normalize qc signal
qc = normFac*qc;

%Compute fitness
ssrVal = -innerprodpsd(params.dataY,qc,params.samplFreq,params.psdVec)^2;

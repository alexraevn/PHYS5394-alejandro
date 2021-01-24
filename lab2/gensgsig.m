function sigVec = genSGSig(dataX,snr,time0,stDev,freq0,phi0)
% Generate a sine-gaussian signal
% S = GENSGSIG(X,SNR,TIME0,STDEV,FREQ0,PHI0)
% Generates a sine-gaussian signal S.
% X is the vector of time stamps at which the value of the signal is to be computed.
% SNR is the matched filtering signal-to-noise ratio of S.
% TIME0 is the time value for the peak of the curve.
% STDEV is the standard deviation for the gaussian curve.
% FREQ0 is the frequency of the sinusoidal oscillation.
% PHI0 is the initial phase of the oscillation.
% S(t) = exp(-(t-time0)^2/(2*stdev^2))*sin(2*pi*freq0*t+phi0)

% Alejandro Reyes, Jan 2021

phaseVec = 2*pi*freq0*dataX+phi0;
gausSig = exp(-(dataX-time0).^2/(2*stDev^2));
sinSig = sin(phaseVec);
sigVec = gausSig.*sinSig;
sigVec = snr*sigVec/norm(sigVec);

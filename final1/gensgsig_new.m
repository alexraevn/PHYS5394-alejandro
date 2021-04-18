function sigVec = gensgsig_new(dataX,snr,P)
% Generate a sine-gaussian signal
% S = GENSGSIG(X,SNR,P)
% Generates a sine-gaussian signal S.
% X is the vector of time stamps at which the value of the signal is to be computed.
% SNR is the matched filtering signal-to-noise ratio of S.
% P.mean is the time value for the peak of the curve.
% P.stdev is the standard deviation for the gaussian curve.
% P.freq0 is the frequency of the sinusoidal oscillation.
% P.phi0 is the initial phase of the oscillation.
% S(t) = exp(-(t-P.mean)^2/(2*P.stdev^2))*sin(2*pi*P.freq0*t+P.phi0)

% Alejandro Reyes, Jan 2021

phaseVec = 2*pi*P.freq0*dataX+P.phi0;
gausSig = exp(-(dataX-P.mean).^2/(2*P.stdev^2));
sinSig = sin(phaseVec);
sigVec = gausSig.*sinSig;
sigVec = snr*sigVec/norm(sigVec);

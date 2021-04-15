%% Test the genSGSig_new function

% Define signal parameter struct
P = struct('mean', 0.5,...
           'stdev', 0.13,...
           'freq0', 22,...
           'phi0', pi/2.0);
snr = 9; % amplitude

% Generate time vector of 1 second
samplFreq = 1024; %hz
samplInterval = 1/samplFreq;
timeVec = 0:samplInterval:1.0; % 1 second in 1024 resolution

% Generate the sine-gaussian signal
sigVec = gensgsig_new(timeVec,snr,P);

% Plot the signal
figure;
plot(timeVec, sigVec, 'Marker', '.', 'MarkerSize', 6);
title('Sine-Gaussian Signal Using Struct');
xlabel('Seconds (s)');
ylabel('Amplitude');

% generate periodogram
nSampl  = length(timeVec); % number of samples
dataLen = timeVec(end)-timeVec(1); % essentially the seconds of data
% kNyq is the frequency sample corresponding to the Nyquist frequency
kNyq = floor(nSampl/2)+1; % Nyquist frequency = nSampl/2*1/dataLen
posFreq = (0:(kNyq-1))*(1/dataLen); % vector of positive Fourier frequencies
% FFT of the signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

% Plot periodogram
figure;
plot(posFreq, abs(fftSig));
title('Periodogram of Sine-Gaussian Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% Plot the sine-gaussian signal
% Signal parameters
time0 = 0.5; % Peak chosen at the center
stDev = 0.15; 
freq0 = 12;
phi0 = pi/2.0;
A = 9;

%% Determination of maximum frequency between sine and gaussian signals
% The frequency function of the gaussian or a(t) function is
% |1/(stDev*sqrt(2*pi)) * exp(t-time0)|
% maxFreq(t) of a(t) occurs at t=time0, time0 ∈ [0,1]
gausMaxFreq = 1/(stDev*sqrt(2*pi));

% The frequency of the sinusoid function is constant
sineMaxFreq = freq0;

% Max frequency is the maximum value between sine and gaussian frequencies
maxFreq = max(gausMaxFreq, sineMaxFreq);

%% Compute sample interval,time vector, and signal
% Sample interval calculation
samplFreq = 5*maxFreq; % sample rate = 5x maximum frequency
samplInterval = 1/samplFreq;

% Time samples vector
timeVec = 0:samplInterval:1.0;
% Number of samples
nSampl = length(timeVec);

% Generate the sine-gaussian signal
sigVec = gensgsig(timeVec,A,time0,stDev,freq0,phi0);

% % Plot the signal
% figure;
% plot(timeVec,sigVec,'Marker','.','MarkerSize',24);
% title('Sine-Gaussian Signal')

%% Plot periodogram
% Length of data
dataLen = timeVec(end)-timeVec(1);
% DFT sample corresponding to Nyquist frequency
kNyq = floor(nSampl/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);
% FFT of signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

% Plot the periodogram
figure;
plot(posFreq, abs(fftSig));
title('Periodogram of Sine-Gaussian Signal')
xlabel('Frequency')
ylabel('Magnitude')


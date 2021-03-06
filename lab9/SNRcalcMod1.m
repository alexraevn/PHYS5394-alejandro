%% How to normalize a signal for a given SNR
% We will normalize a signal such that the Likelihood ratio (LR) test for it has
% a given signal-to-noise ratio (SNR) in noise with a given Power Spectral 
% Density (PSD). [We often shorten this statement to say: "Normalize the
% signal to have a given SNR." ]

% Include functions folder
% addpath ./FUNCTIONS
% SDM **************************
% Why didn't you addpath to DATASCIENCE_COURSE/DETEST? 
% As instructed in the lab (slide #2), the codes were provided. 
addpath ../lab2
%*******************************

% This is the target SNR for the LR
snr = 10;

%% Generate Sine Gaussian Signal
% Data generation parameters
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

% Sine Gaussian signal parameters
time0 = 0.5; % Peak chosen at the center
stDev = 0.08; 
freq0 = 30;
phi0 = pi/2.0;
% Amplitude value does not matter as it will be changed in the normalization
A = 1; 
sigVec = gensgsig(timeVec,A,time0,stDev,freq0,phi0);

%% Use noisePSD inline function to generate the PSD vector 
% We will use the noise PSD used in colGaussNoiseDemo.m but add a constant
% to remove the parts that are zero. (Exercise: Prove that if the noise PSD
% is zero at some frequencies but the signal added to the noise is not,
% then one can create a detection statistic with infinite SNR.)
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;

% Generate the PSD vector to be used in the normalization. Should be
% generated for all positive DFT frequencies. 
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);
% figure;
% plot(posFreq,psdPosFreq);
% axis([0,posFreq(end),0,max(psdPosFreq)]);
% xlabel('Frequency (Hz)');
% ylabel('PSD ((data unit)^2/Hz)');

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,psdPosFreq);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Test
%Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,psdPosFreq);
end
%Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    % Add normalized signal
    dataVec = noiseVec + sigVec;
    llrH1(lp) = innerprodpsd(dataVec,sigVec,sampFreq,psdPosFreq);
end
%%
%Signal to noise ratio estimate
estSNR = (mean(llrH1)-mean(llrH0))/std(llrH0);

figure;
histogram(llrH0);
hold on;
histogram(llrH1);
xlabel('LLR');
ylabel('Counts');
legend('H_0','H_1');
title(['Estimated SNR = ',num2str(estSNR)]);

%% Plot data and signal realizations
figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sigVec);
title('Data and Signal Realizations');
xlabel('Time (sec)');
ylabel('Data');

%% Generate and plot the periodogram of signal and data

% FFT of the data
fftData = fft(dataVec);
% Take positive fourier frequencies
fftData = fftData(1:kNyq);
% FFT of signal (positive frequencies)
fftSig = fft(sigVec);
fftSig = fftSig(1:kNyq);

% Plot periodogram of signal
figure;
plot(posFreq, abs(fftData));
hold on;
plot(posFreq, abs(fftSig));
title('Periodogram of Signal and Data');
xlabel('Frequency (Hz)');
ylabel('Periodogram');

%% Generate and plot the spectrogram

% Set the window length and overlap
winLen = 0.09; %s
ovrlp  = 0.07; %s

% Convert to integer number of samples
winLenSampl = floor(winLen*sampFreq);
ovrlpSampl  = floor(ovrlp*sampFreq);
[S,F,T]     = spectrogram(dataVec, winLenSampl, ovrlpSampl, [], sampFreq);

% Plot spectrogram
figure;
imagesc(T, F, abs(S)); axis xy;
title('Spectrogram of Data');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

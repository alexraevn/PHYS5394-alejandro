%% How to normalize a signal for a given SNR
% We will normalize a signal such that the Likelihood ratio (LR) test for it has
% a given signal-to-noise ratio (SNR) in noise with a given Power Spectral 
% Density (PSD). [We often shorten this statement to say: "Normalize the
% signal to have a given SNR." ]

% Include functions folder
%addpath ./FUNCTIONS
% SDM **************************
% Why didn't you addpath to DATASCIENCE_COURSE/DETEST? 
% As instructed in the lab (slide #2), the codes were provided. 
addpath ../lab2
%*******************************
addpath ./DATA

% This is the target SNR for the LR
snr = 10;

% Data generation parameters
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%% Generate the signal that is to be normalized
% Sine Gaussian signal parameters
time0 = 0.5; % Peak chosen at the center
stDev = 0.08; 
%SDM ******************************
%We want to look at signals in the sensitive band [50, 700] Hz
freq0 = 200;
%***********************************
%freq0 = 20;
phi0 = pi/2.0;
% Amplitude value does not matter as it will be changed in the normalization
A = 1; 
sigVec = gensgsig(timeVec,A,time0,stDev,freq0,phi0);

%% Use iLIGOSensitivity.m to generate a PSD vector
% (Exercise: Prove that if the noise PSD is zero at some
% frequencies but the signal added to the noise is not,
% then one can create a detection statistic with infinite SNR.)

% Read the LIGO sensitivity file
LIGOSen = load('iLIGOSensitivity.txt', '-ascii');

% Manual pre-filtering: set S(f<=lowCutoff)=S(f=lowCutoff)
lowCutoff  = 50; %Hz
% LIGOSen(42,1) = 49.79 Hz, we'll use index 42 for the low cutoff of 50 Hz
lowSen = LIGOSen(42,2); %S value around f = lowCutoff
LIGOSen(1:41,2) = lowSen;

% End vector at highcutoff = sampFreq/2
% LIGOSen(67:1) = 490.37 Hz, index 67 for high cutoff of 512 Hz
% This is al alteration from 490 to 512 Hz in the sensitivity
% I tested polyfit to get a value at exactly 512 Hz but have been unsuccessful
highCutoff = sampFreq/2; %Hz
LIGOSen(67,1) = highCutoff;
LIGOSen = LIGOSen(1:67,:); % Truncate

% Insert f = 0 sampling freq
LIGOSen(2:end+1,:) = LIGOSen;
LIGOSen(1,1) = 0; %Hz
LIGOSen(1,2) = lowSen;

% Interpolate sensitivity values to match length(psdVals) = kNyq = floor(nSamples/2)+1
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen); % Positive DFT frequencies

psdVals = interp1(LIGOSen(:,1),LIGOSen(:,2),posFreq);

%SDM****************************
%Both statgaussnoisegen and innerprodpsd (in DATASCIENCE_COURSE) use PSD,
%not sqrt(PSD), which is what is provided in iLIGOSensitivity.txt). Review
%the lab where you simulated LIGO noise.
psdVals = psdVals.^2;
%*******************************

% figure;
% loglog(posFreq,psdVals);
% axis([0,LIGOSen(1:end,1),0,LIGOSen(1:end,2)]);
% xlabel('Frequency (Hz)');
% ylabel('PSD ((data unit)^2/Hz)');

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,psdVals);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Test
%Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:), psdVals(:)],100,sampFreq);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,psdVals);
end
%Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:), psdVals(:)],100,sampFreq);
    % Add normalized signal
    dataVec = noiseVec + sigVec;
    llrH1(lp) = innerprodpsd(dataVec,sigVec,sampFreq,psdVals);
end
%%
% Signal to noise ratio estimate
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




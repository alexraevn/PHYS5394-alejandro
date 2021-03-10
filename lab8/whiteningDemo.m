%% Read testData, whiten signal and plot spectrogram before and after
% Read testData.txt
data = load('testData.txt');
% Split file into times and data vectors
sampTimes = data(:,1);
dataVals  = data(:,2);
% Calculate sample frequency
sampFreq = length(sampTimes)/round(sampTimes(end)); % This works for testData.txt as it is an integer number of seconds long
% First 5 seconds of data
colGausTimes = sampTimes(1:5*sampFreq);
colGaus      = dataVals(1:5*sampFreq);

%% Estimate Power Spectral Density function
% estimate PSD of first 5 seconds of data using pwelch, skipping first 0.25 seconds of silence
[pxx,f] = pwelch(colGaus, 80,[],[],sampFreq);
% colGaus is first 5 seconds of data in 'testData.txt'
% 80 is the length of the Hamming window which gives a smooth PSD (5120/80 = 64 windows)
% 50% window overlap and 256 FFT points set by default

%% Design Whitening Filter
% FIR2 filter for target T(f)=sqrt(PSD) where PSD is the estimated PSD of first 5 seconds
fltrOrdr = 500;
b = fir2(fltrOrdr,f/(sampFreq/2),sqrt(pxx));
% Whiten testData
whitenedData = sqrt(sampFreq)*fftfilt(b,dataVals); % Normalized to two sided PSD

%% Spectrograms
% Declare length of windows and their overlap
winLen = 0.08;
ovrlp  = 0.07;

% Convert to integer number of samples
winLenSamp = floor(winLen*sampFreq);
ovrlpSamp  = floor(ovrlp*sampFreq);

% Generate spectrogram of testData.txt
[Sin,Fin,Tin] = spectrogram(dataVals, winLenSamp, ovrlpSamp, [], sampFreq);

% Plot data spectrogram
figure;
imagesc(Tin, Fin, abs(Sin)); axis xy;
title('Spectrogram of testData');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Generate spectrogram of whitened data
[Sout,Fout,Tout] = spectrogram(whitenedData, winLenSamp, ovrlpSamp, [], sampFreq);

% Plot whitened data spectrogram
figure;
imagesc(Tout, Fout, abs(Sout)); axis xy;
title('Spectrogram of Whitened testData');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

figure;
subplot(2,1,1), plot(sampTimes, dataVals), title('testData');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2), plot(sampTimes, whitenedData), title('Whitened Data');
xlabel('Time (s)');
ylabel('Amplitude');

% Testing alignment truncating first 0.25 s of silence
% figure;
% subplot(3,1,1), plot(sampTimes, dataVals);
% subplot(3,1,2), plot(sampTimes(1:end-250), whitenedData(251:end));
% subplot(3,1,3), plot(sampTimes(1:end-250), whitenedData(251:end)-dataVals(1:end-250));


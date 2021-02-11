%% Generate 3 sinusoids
% Define signal parameters
A1 = 10; % SNR
f1 = 100; % Sine frequency
phi1 = 0; % Initial phase

A2 = 5;
f2 = 200;
phi2 = pi/6;

A3 = 2.5;
f3 = 300;
phi3 = pi/4;

% Generate time vector
samplFreq = 1024; % Hz
maxFreq = 0.5*samplFreq % Same as 1.0 = half the sample rate required by fir1
nSamples = 2*samplFreq;
timeVec = (0:(nSamples-1))/samplFreq;

% Declare the 3 sinusoids
s1 = A1*sin(2*pi*f1*timeVec+phi1);
s2 = A2*sin(2*pi*f2*timeVec+phi2);
s3 = A3*sin(2*pi*f3*timeVec+phi3);

% Add 3 sinusoids
s = s1 + s2 + s3;

%% Filter design
% Design filter 1 to let s1 pass (lowpass)
filtOrder = 30; % Use same order for all
lPass = fir1(filtOrder, 150/maxFreq, 'low'); % Pass frequencies below 150 Hz
% Apply the filter
filtS1 = fftfilt(lPass, s);

% Design filter 2 to let s2 pass (bandpass)
hBandPass = fir1(filtOrder, 150/maxFreq, 'high'); % Pass frequencies above 150 Hz
lBandPass = fir1(filtOrder, 250/maxFreq, 'low'); % Pass frequencies below 250 Hz
% Apply bandpass
filtS2 = fftfilt(hBandPass, s); % low cutoff
filtS2 = fftfilt(lBandPass, filtS2); % high cutoff

% Design filter 3 to let s3 pass (highpass)
hPass = fir1(filtOrder, 250/maxFreq, 'high'); % Pass frequencies above 250 Hz
% Apply highpass
filtS3 = fftfilt(hPass, s);

%% Periodograms
% Length of data
dataLen = timeVec(end)-timeVec(1);
% DFT sample corresponding to Nyquist frequency
kNyq = floor(length(timeVec)/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);

% FFT of unfiltered signal s
fftSig = fft(s);
fftSig = fftSig(1:kNyq);

% FFT of filtered signal 1
fftSig1 = fft(filtS1);
fftSig1 = fftSig1(1:kNyq);

% FFT of filtered signal 2
fftSig2 = fft(filtS2);
fftSig2 = fftSig2(1:kNyq);

% FFT of filtered signal 3
fftSig3 = fft(filtS3);
fftSig3 = fftSig3(1:kNyq);

%% Plot the periodograms
fig = figure;

subplot(2,3,[1,3]); % Plot the periodogram of the input signal once
hold on;
plot(posFreq, abs(fftSig));
title('Input signal periodogram');

subplot(2,3,4); % Plot a periodogram of each filtered signal
hold on;
plot(posFreq, abs(fftSig1));
title('Lowpassed output periodogram');

subplot(2,3,5);
hold on;
plot(posFreq, abs(fftSig2));
title('Bandpassed output periodogram');

subplot(2,3,6);
hold on;
plot(posFreq, abs(fftSig3));
title('Highpassed output periodogram');

hand = axes(fig, 'visible', 'off'); % All of this just to set up common axis labels
hand.Title.Visible='on';
hand.XLabel.Visible='on';
hand.YLabel.Visible='on';
ylabel(hand, 'Magnitude');
xlabel(hand, 'Frequency (Hz)');

%% Plot the signals
% Extra plots of filtered signals imposed over input
fig2 = figure;

subplot(3,1,1); % 3 row, 1 column figure of 3 plots
hold on;
plot(timeVec, s);
plot(timeVec, filtS1);
title('Lowpass cutoff 150 Hz (s1)');

subplot(3,1,2);
hold on;
plot(timeVec, s);
plot(timeVec,filtS2);
title('Bandpass cutoff 150-250 Hz (s2)');

subplot(3,1,3);
hold on;
plot(timeVec, s);
plot(timeVec, filtS3);
title('Highpass cutoff 250 Hz (s3)');

hand = axes(fig, 'visible', 'off');
hand.Title.Visible='on';
hand.XLabel.Visible='on';
hand.YLabel.Visible='on';
ylabel(hand, 'Magnitude');
xlabel(hand, 'Frequency (Hz)');

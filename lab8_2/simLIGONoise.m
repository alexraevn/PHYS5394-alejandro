%% Simulate LIGO Noise
% Read the initial LIGO sensitivity PSD
LIGOSen = load('iLIGOSensitivity.txt', '-ascii');

%% Manual pre-filtering: set S(f<=50)=S(f=50) and S(f>=700)=S(f=700)

lowCutoff  = 50; %Hz
% lowSen     = 
highCutoff = 700; %Hz
% highSen    =

% LIGOSen(42,1) = 49.79 Hz, we'll use index 42 for the low cutoff of 50 Hz
lowSen = LIGOSen(42,2); %S value around f = lowCutoff

for i = 1:41
    LIGOSen(i,2) = lowSen;
end

% LIGOSen(70:1) = 717.95 Hz, index 70 for high cutoff of 700 Hz
highSen = LIGOSen(70,2);

for i = 71:length(LIGOSen(:,1))
    LIGOSen(i,2) = highSen;
end

%% Insert f = 0, f = 1/2 sampling freq
% f=0
LIGOSen(2:end+1,:) = LIGOSen;
LIGOSen(1,1) = 0; %Hz
LIGOSen(1,2) = lowSen;

% Since we don't need frequencies above 700 Hz, I will truncate at 3 kHz
% LIGOSen(85,1) = 2999.1 Hz, so stop at index 85 and round to 3 kHz
% Set the sampling frequency at 6 kHz for simplicity
sampFreq = 6000;

LIGOSen = LIGOSen(1:85,:); % Truncate
LIGOSen(85,1) = 3000; %Hz

% loglog(LIGOSen(:,1),LIGOSen(:,2))

%% Generate WGNoise realization and pass through LIGO PSD filter
% Number of samples = sampFreq * 'time' length. Using 2 seconds:
nSamples = 2*sampFreq;
% Generate the Gaussian noise vector
noise = randn(1,nSamples);

% Design FIR filter
freqVals  = LIGOSen(:,1);
psdVals   = LIGOSen(:,2);
filtrOrdr = 600; % Using 600 order and a similar pwelch window consistently reduces the jaggedness
                 % at the frequencies of interest.

b = fir2(fltrOrdr,freqVals/(sampFreq/2),sqrt(psdVals));

noiseLIGOReal = sqrt(sampFreq)*fftfilt(b,noise);

%% Estimate the PSD using pwelch

[pxx, f] = pwelch(noiseLIGOReal, 600, [], [], sampFreq);

figure;
plot(f/1000,pxx);
title('Estimated PSD (pwelch) of LIGO Noise Realization');
xlabel('Frequency (kHz)');
ylabel('PSD');

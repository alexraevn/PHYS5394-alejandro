%% How to normalize a signal for a given SNR
% We will normalize a signal such that the Likelihood ratio (LR) test for it has
% a given signal-to-noise ratio (SNR) in noise with a given Power Spectral 
% Density (PSD). [We often shorten this statement to say: "Normalize the
% signal to have a given SNR." ]

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
freq0 = 20;
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

% Generate the PSD vector to be used in the normalization. Should be
% generated for all positive DFT frequencies.
% psdFreq = LIGOSen(:,1);
% psdVals = LIGOSen(:,2);

figure;
loglog(posFreq,psdVals);
% axis([0,LIGOSen(1:end,1),0,LIGOSen(1:end,2)]);
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

%%
% A noise realization
figure;
plot(timeVec,noiseVec);
xlabel('Time (sec)');
ylabel('Noise');
%%
% A data realization
figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sigVec);
xlabel('Time (sec)');
ylabel('Data');




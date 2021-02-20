%% Generate and plot the sine-gaussian signal
% Signal parameters
time0 = 0.5; % Peak chosen at the center
stDev = 0.15; 
freq0 = 9;
phi0 = pi/2.0;
A = 9;

% The maximum frequency is equal to the constant frequency of the sinusoid
maxFreq = freq0;

% Package signal parameters into a cell array to feed into genPlotSig function
argsVector = {time0, stDev, freq0, phi0, A, maxFreq};

%% Generate and plot signals
% Signal sampled at 5*maximum frequency
genPlotSig(5,argsVector)

% Signal smapled at 0.5*maximum frequency
genPlotSig(0.5,argsVector)

%% Function to generate time vector, signal, and plot using different multiples of maximum frequency
function genPlotSig(samplFactor,argsVec)
% Generate signal for a given factor of maximum frequencies and plot it
% SAMPLFACTOR is the multiple of maximum frequencies signal is sampled at
% argsVec is a cell array of arguments in the following order
%   argsVec = {time0, stDev, freq0, phi0, A, maxFreq}

% Compute sample frequency, interval, and time vector
samplFreq = samplFactor*argsVec{6};
samplInterval = 1/samplFreq;
timeVec = 0:samplInterval:1.0;

% Generate the sine-gaussian signal
sigVec = gensgsig(timeVec,argsVec{5},argsVec{1},argsVec{2},argsVec{3},argsVec{4});

% Plot the signal
figure;
plot(timeVec,sigVec,'Marker','.','MarkerSize',24);
title(sprintf('Sine-Gaussian Signal, Sampling at %.1gx Maximum Frequency',samplFactor))
xlabel('Magnitude')
ylabel('Time sample')
end
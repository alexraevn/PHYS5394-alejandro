%% Generate dataX vector, P struct and pass to gensgsig_new to plot for SNR 10, 12, 15.

% Generate time series
samplFreq = 1024; % Hz
% dataLen = 4.0; % seconds
samplInterval = 1/samplFreq;
dataX = 0:samplInterval:1.0; % time series vector

% Generate the struct of sine-gaussian signal parameters
P = struct('mean', 0.5,...
           'stdev', 0.13,...
           'freq0', 22,...
           'phi0', pi/2.0);

% Function handle to call genSGSig_new
H = @(snr) gensgsig_new(dataX, snr, P);

%% Plot sine-gaussian signal at 10, 12, and 15 SNR
% Specify the needed snr(s) in scalable array S
S = [10 12 15];

% Generate legend cell array and new figure
legVec = [];
figure;
% Loop through snrs, generate signal and add to plot
for snr = S
    % Append snr to legend cell array
    legVec{end+1} = ['snr=',num2str(snr)];
    sigVec = H(snr);
    % Plot the signal
    plot(dataX, sigVec);
    hold on;
end
% Display title, axis labels, and legend
title('Sine-Gaussian Signal at Various SNRs');
xlabel('Time (s)');
ylabel('Amplitude');

hold off;
legend(legVec);

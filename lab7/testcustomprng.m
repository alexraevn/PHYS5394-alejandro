%% Test the customrand and customrandn functions
%% Generate probability distributions
% Set the number of trials for each probability distribution
trials = 10000;

for i = 1:trials;
% Generate the uniform PDF vector of length TRIALS
uPDF(i) = customrand(-2,1); % Between -2,1

% Generate the normal PDF vector of length TRIALS
nPDF(i) = customrandn(1.5,2.0); % Mean = 1.5, Std = 2.0

end

%% Plot PDFs with histograms
% The histograms are normalized by probability density function estimate
figure;
subplot(2,1,1), plot(uPDF), title('Uniform PDF'), xlabel('Trial'), ylabel('Value');
subplot(2,1,2), histogram(uPDF, 'Normalization', 'pdf'), title('Histogram'), xlabel('Value'), ylabel('Probability');
figure;
subplot(2,1,1), plot(nPDF), title('Normal PDF'), xlabel('Trial'), ylabel('Value');
subplot(2,1,2), histogram(nPDF, 'Normalization', 'pdf'), title('Histogram'), xlabel('Value'), ylabel('Probability');

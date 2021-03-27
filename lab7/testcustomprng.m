%% Test the customrand and customrandn functions
%% Generate probability distributions
% Set the number of trials for each trial run
trials = 10000;
% Declare PDF parameters
a = -2;
b = 1;
m = 1.5; %mean
s = 2.0; %std

for i = 1:trials;
% Generate the vector of trials from uniform PDF
uniformTrials(i) = customrand(a,b); % Between -2,1

% Generate the vector of trials from normal PDF
normalTrials(i) = customrandn(m,s); % Mean = 1.5, Std = 2.0

end

%% Plot PDFs with histograms
% Uniform PDF
uLimits = [a-1,b+1]; % Set x axis limits for plots
syms x1;
uPDF = piecewise(uLimits(1)<=x1<-2, 0, a<=x1<=b, 1/(b-a), 1<x1<=uLimits(2), 0); % = 1/(b-a), [-2,1], = 0, otherwise
% Normal PDF
nLimits = [m-4*s,m+4*s]; % nPDF(x2 = m+-4*s) <= 0.035, which looks good in most cases 
syms x2;
nPDF = (1/[sqrt(2*pi)*s])*exp(-[x2-m]^2/[2*s^2]);

% Two figures with 2 subplots each
figure;
hold on;
histogram(uniformTrials, 'Normalization', 'pdf'), xlim(uLimits), title(sprintf('Histogram of %d Trials', trials)), xlabel('Value'), ylabel('Probability');
fplot(uPDF), title('Uniform PDF, limits(-2,1)'), xlabel('Trial'), ylabel('Value');
figure;
hold on;
histogram(normalTrials, 'Normalization', 'pdf'), xlim(nLimits), title(sprintf('Histogram of %d Trials',trials)), xlabel('Value'), ylabel('Probability');
fplot(nPDF, nLimits), title('Normal PDF, Mean=1.5, Std=2.0'), xlabel('Trial'), ylabel('Value');


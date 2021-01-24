%% Plot the sine-gaussian signal
% Signal parameters
time0 = 0.5;
stDev = 0.15;
freq0 = 12;
phi0 = pi/2.0;
A = 9;

% Sample interval calculation
samplFreq = 7*sqrt(5)*freq0; % Irrational number used to remove cyclic patterns in datapoints
samplInterval = 1/samplFreq;

% Time samples vector
timeVec = 0:samplInterval:1.0;

% Generate the sine-gaussian signal
sigVec = gensgsig(timeVec,A,time0,stDev,freq0,phi0);

% Plot the signal
figure;
plot(timeVec,sigVec,'Marker','.','MarkerSize',24);


%% Test harness for CRCBPSO 
% The fitness function called is CRCBPSO_GRIEWANK. 
ffparams = struct('rmin',-600,...
                     'rmax',600 ...
                  );
% Fitness function handle.
fitFuncHandle = @(x) crcbpso_griewank(x,ffparams);
%%
% Define structure P for customizing crcbpso
P = struct('popSize',[],...
           'maxSteps',[],...
           'boundaryCond',[],...
           'nbrhdSz',[]...
          );
% Number of dimensions
nDim = 30;

% Reset rng seed
rng('default')


%% Set up multiple iterations

% Number of pso trials
trials = 4;
bestLocations = zeros(trials,nDim); % Vector of best standard coordinates
bestFitnesses = zeros(trials,1); % Vector of best fitness values

% parpool(4);

for n = 1:trials
    % Call crcbpso
    psoOut = crcbpso(fitFuncHandle,nDim,P);
    % Append coord and fitness from psoOut
    bestLocations(n,:) = psoOut.bestLocation;
    bestFitnesses(n) = psoOut.bestFitness;
end

% Ultimate fitness value and its index
[ultFitness,ultIndex] = min(bestFitnesses);
% Transform standard coordinate to real coordinate
% Ultimate location in real coordinates
[~,ultLocation] = fitFuncHandle(bestLocations(ultIndex,:));

%% Display
disp(['Best location:',num2str(ultLocation)]);
disp(['Best fitness:', num2str(ultFitness)]);
%% Estimated parameters
% Best standardized and real coordinates found.
% stdCoord = psoOut.bestLocation;
% [~,realCoord] = fitFuncHandle(stdCoord);


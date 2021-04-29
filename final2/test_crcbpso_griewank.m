%% Test harness for CRCBPSO 
% The fitness function called is CRCBPSO_GRIEWANK. 
ffparams = struct('rmin',-600,...
                     'rmax',600 ...
                  );
% Fitness function handle.
fitFuncHandle = @(x) crcbpso_griewank(x,ffparams);
%%
% Define structure P for customizing crcbpso
P = struct('popSize',40,... % Number of particles
           'maxSteps',4000,... % Number of steps
           'boundaryCond',[],...
           'nbrhdSz',[]...
          );
% Number of dimensions
nDim = 30;

% Reset rng seed
% rng('default')

%% Run multiple crcbpso

% Number of pso trials
trials = 6; % Low number of trials which yields precise results
% Define empty arrays of locations / fitnesses
bestLocations = zeros(trials,nDim); % Vector of best standard coordinates
bestFitnesses = zeros(trials,1); % Vector of best fitness values

% Use Parallel Computing Toolbox parfor loop
% Set up number of workers using parloop(n)
parfor n = 1:trials
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

%% Run information
% Average fitness found in all trials
avgFitness = mean(bestFitnesses);
% Standard deviation
stdFitness = std(bestFitnesses);

% particles = P.popSize
steps = P.maxSteps

%% Display
disp('%========%')
% disp([num2str(trials), ' trials, 40 particles, 2000 steps']) % Trial run
% disp([num2str(particles), ' particles, 72 trials, 2000 steps']) % Particle run
disp([num2str(steps), ' steps,', num2str(P.popSize),' particles,',num2str(trials),' trials'])
disp(['Best fitness:', num2str(ultFitness)]);
disp(['Average fitness:', num2str(avgFitness)]);
disp(['Standard deviation:', num2str(stdFitness)]);
disp('%========%')

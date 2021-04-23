%% Plot benchmark fn
% Set x,y limits and step resolution
lim = 50;
res = 0.1;

% Function handles to test
fn_ackley = @(x,n) -20*exp(-0.2*sqrt(x.^2/n)) - exp(cos(2*pi*x)) + 20 + exp(1);

fn_griewank = @(x,n) x.^2/4000 - cos(x./sqrt(n)) + 1;

% Define X, Y, Z
[X,Y] = meshgrid(-lim:res:lim,-lim:res:lim);

% Z = sqrt(X.^2 + Y.^2); %Sphere
% Z = fn_ackley(X,1) + fn_ackley(Y,2);
Z = fn_griewank(X,1) + fn_griewank(Y,2);
% Z = 100*(Y-X.^2).^2 + (1-X).^2; %Rosenbrock % It does not look nice


surf(X,Y,Z,'FaceAlpha', 0.6,'EdgeColor','none');
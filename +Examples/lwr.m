function data = lwr( )

global data
data.Vvec = [];
data.Hvec = [];
data.u_a = [];
data.u_b = [];
data.delta = [];

%% Set configurations
conf = Configuration();

conf.model = Model.LWR;
conf.solver = Flux.Rusanov;

conf.timeInt = @TimeIntegration.FE;
conf.tMax = 30;
conf.CFL = 0.1;

conf.bc = Mesh.BC.ControlledStabilityAndInvarianceSeparate;

conf.mesh = Mesh.Cartesian([-1,1], 100);

conf.initial = @(x) 0.1 - 0.1 * sin(pi * x);


%% Run solver
soln = runSolver(conf);


%% Display data
Plot.plotSolution(soln, [0 1]);

data.soln = soln;

end
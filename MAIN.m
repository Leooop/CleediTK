
% Except for the temperature that must be given in °C, every other dimensional parameter must be supplied in its international units.

%% OUTPUTS PARAMETERS %%
OPT.write_PT_path = 1; % write P-T path in file
OPT.color_phase = 9; % The phase to associate with lines color on plots
OPT.hold.fig2 = 0; % hold figure 2 [yes=1/no=0]
OPT.hold.fig3 = 0; % same as above
OPT.hold.fig4 = 0;

%% DYNAMICS PARAMETERS %%
year = 3600*24*365.25; % a year in seconds
DYN.v_slab = 0.005/year; % (m/s) slab velocity
DYN.dt = 3e11; % (s) timestep
               % be careful to satisfy the following for numerical stability :  dt < min(min_caracteristic_reaction_time) 
               % caracteristic_reaction_time = 1/min(reaction_constant)
               % no need to constrain on thermal diffusivity because T
               % field is at equilibrium.

%% KINETIC PARAMETERS %%

KIN.E = 1e5; % activation energy [J.mol-1]
KIN.R = 8.14; % perfect gas constant
KIN.phases_diagram = 'facies.mat'; % matlab file containing the phase diagram along with associated T and P vectors.
[KIN.get_phase, KIN.nphases, KIN.Pdiag, KIN.Tdiag, KIN.phases_diag] = load_phase_diagram(KIN.phases_diagram); % get_phase(T,P) is an callable interpolation object returning the phase index at (T,P) and nphases is the total number of phases

% reaction constants matrix, 
% the component K_ij contains the constant associated with the reaction from phase i to phase j.
% Its size must be nphases^2 and diagonal terms should be zero.

% example prescribing a similar eclogitization reaction constant from all facies, and setting every other
% constant to one tenth of eclogitization constant. Symmetric matrix.
nphases = KIN.nphases;
k_eclogite = 2.27e-10; % fit baxter for eclogitization
k_mat = ones(nphases,nphases) .* k_eclogite .* 10; 
k_mat(:,9) = k_mat(:,9) ./ 10; % eclogite is phase 9. column 9 contains eclogitization constants from all phases.
for i = 1:nphases
    k_mat(i,i) = 0; % zero diagonal
end
KIN.k_mat = k_mat;

%% GEOMETRIC PARAMETERS %%
GEOM.x0    = [-100e3 440e3]; %width of box, 0 = trench location
GEOM.y0    = [0 150e3]; %length of box
GEOM.res   = [100,50]; %resolution of the map in x,y
GEOM.dx = diff(GEOM.x0)/GEOM.res(1);
GEOM.dy = diff(GEOM.y0)/GEOM.res(2);
GEOM.slab_dip = 20; % degree
%GEOM.H_slab = 0;

GEOM.T_dependent_lplate_thickness = 0 ; % [0/1] to use either 'lplate_thickness' or 'T_lplate_thickness' to determine rigid subducting plate thickness.
GEOM.T_lplate_thickness = 0 ; %(°C) temperature at the bottom of the subducting lithosphere, used to determine a constant thickness from the left boundary temperature profile.
GEOM.lplate_thickness = 0 ; % (m) static upper plate thickness.
GEOM.uplate_thickness = 30e3 ; % (m) static upper plate thickness, zero velocity.
GEOM.radiogenic_uplate_thickness = 30e3; % (m) thickness of the radiogenic upper plate
GEOM.radiogenic_lplate_thickness = 20e3; % (m) thickness of the radiogenic lower (subducting) plate

%% THERMAL PARAMETERS %%
THERMAL.heat_prod = 3e-13; % [K/s] : puissance volumique [W/m^3] divisée par rho*Cp [kg*m^(-3)*J*kg^(-1)*K^(-1)]
THERMAL.diffusivity = 1e-6; % (m^2/s) thermal diffusivity

%% BOUNDARY CONDITIONS %%
% These BC definitions are not particularly consistent and should be
% redesigned for clarity

% boundary types : 
%   - 1 : T imposed constant along the boundary
%   - 3 : mixte (dirichlet and Neumann) (T imposed where inflow, grad(T)=0 orthogonaly to the boundary where outflow => no diffusive outflow)
BC.top.type = 1; % constant T
BC.bot.type = 3; % no diffusive outflow
BC.left.type = 3; % half space cooling geotherm using BC.age_slab
BC.right.type = 3; % linear vertical T gradient from BC.top.value to BC.right.value in the inflow part, no diffusive outflow deeper.

BC.top.val = 0; % T value top, also used to constrain minimum value at the left boundary
BC.bot.val = 1300; % T inflow value bottom, also used to constrain maximum value at the left boundary
BC.left.val = 0; % not used...
BC.right.val = 1300; % T at the limit between inflow and outflow.

BC.age_slab = 80e6*year ;


%% RUN CODE %%
%%%%%%%%%%%%%%

[P,T,X_mat] = thermokin(OPT,DYN,GEOM,KIN,THERMAL,BC);


%
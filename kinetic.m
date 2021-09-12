function [X_mat_RK, id_vec, t_vec] = kinetic(P,T,K_PARAMS,DYN)
%% get parameters 
dt = DYN.dt;
E = K_PARAMS.E; 
R = K_PARAMS.R; 
n_id = K_PARAMS.nphases;
get_phase = K_PARAMS.get_phase;
k_mat = K_PARAMS.k_mat;

% initial conditions
X0 = zeros(1,n_id);
X0(get_phase(T(1),P(1))) = 1; % composition corresponds to our location in the facies 

% ponderation matrix
pond_mat = ones(n_id,n_id); % pondération de concentrations 
for i = 1:n_id 
    pond_mat(i,i) = 0; % zero diagonal
end

% constrain dt according to maximum reaction rate :

%% integration routine

niter = length(T); % nb d'iteration
X_mat_RK = zeros(niter,n_id);
init_phase = get_phase(T(1),P(1));
X_mat_RK(1,init_phase) = 1; 
id_vec = zeros(1,niter);
id_vec(1) = get_phase(T(1),P(1));
t_vec = 0:dt:dt*(niter-1);

for it = 2:niter
    current_phase = get_phase(T(it),P(it));
    id_vec(it) = current_phase;
    for pid = 1:n_id
        %X_mat(it,pid) = X_mat(it-1,pid) * (1 - k_mat(pid,current_phase).*exp(-E/(R*(T(it)+ 273.15)))*dt);
        X_mat_RK(it,pid) = RK4_step(X_mat_RK(it-1,pid),k_mat(pid,current_phase),T(it)+273.15,E,R,dt); % RK4
    end
    %X_mat(it,current_phase) = 1 - pond_mat(current_phase,:)*X_mat(it,:)';
    X_mat_RK(it,current_phase) = 1 - pond_mat(current_phase,:)*X_mat_RK(it,:)';
end

end
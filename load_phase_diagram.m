function [get_phase, n_id, Pmat, Tmat, phases_mat] = load_phase_diagram(file)

load(file)
phases_mat = double(id_facies);
[Pmat,Tmat] = ndgrid(P_facies,T_facies);
Pup = flipud(Pmat);
idup = flipud(phases_mat);
get_phase = griddedInterpolant(Tmat',Pup',idup','nearest');
%get_phase_metam = @(P,T) F(T,P);

n_id = length(unique(phases_mat));

end
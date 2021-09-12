function u_next = RK4_step(u,k,T,E,R,dt)

% Arhenius law :
du = @(u) - k * exp(-E/(R*T)) * u;

k1 = du(u);
k2 = du(u + dt*k1/2);
k3 = du(u + dt*k2/2);
k4 = du(u + dt*k3);

u_next = u + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
end
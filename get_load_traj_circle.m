 function load_traj = get_load_traj_circle(t)
%% Circular Displacement  

r = 3; % radius
w = 1; % angular velocity
c = 0; % center

%% Desired Load Position Generation
load_traj.xL = c + r*[sin(w*t); cos(w*t); 1];
load_traj.dxL = r*[w*cos(w*t); -w*sin(w*t); 0];
load_traj.d2xL = r*[-w^2*sin(w*t); -w^2*cos(w*t); 0];
load_traj.d3xL = r*[-w^3*cos(w*t); w^3*sin(w*t); 0];
load_traj.d4xL = r*[w^4*sin(w*t); w^4*cos(w*t); 0];
load_traj.d5xL = r*[w^5*cos(w*t); -w^5*sin(w*t); 0];
load_traj.d6xL = r*[-w^6*sin(w*t); w^6*cos(w*t); 0];

end


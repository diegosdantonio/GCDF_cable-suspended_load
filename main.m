% Geometric Tracking Control of a Quadrotor UA V on SE(3) , 
% Taeyoung Lee, Melvin Leok, and N. Harris McClamroch 

% Date: June-12-2020
% Last Updated: June-12-2020

close all;

addpath('./Geometry-Toolbox');

%% Parameters
data.params.mQ = 0.5;
data.params.mL = 0.1;
data.params.J = [2.32e-3,0,0;0,2.32e-3,0;0,0,4e-3];
data.params.g = 9.81;
data.params.e1 = [1;0;0];
data.params.e2 = [0;1;0];
data.params.e3 = [0;0;1];
data.params.l = 1;

%% initial condition

xL = [-3;-3;2];
vL = zeros(3,1);
th = 125*pi/180;
q = [-sin(th);0;cos(th)];
omega = [0;0;0];
R = eye(3,3);
Omega = [0;0;0];


x_0 = [xL; q; reshape(R, 9,1); vL ;omega; Omega];


%% ODE Dynamical Equations solver

odeopts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t, x] = ode15s(@odefun_control, [0 10], x_0, odeopts, data);


%% Compute  quantities

disp('Computing State Variables and Configuration Errors ') ;
ind = round(linspace(1, length(t), round(0.1*length(t)))) ;
for j=ind
    [~,xLd_, qd_, Rd, vL_ ,omega_, Omega_,  f_, M_] = odefun_control(t(j), x(j,:)', data) ;
    xLd(j,:)=xLd_';
    M(j,:) = M_';
    f(j,:) = f_';
end

%% Plot
ind = round(linspace(1, length(t), round(0.1*length(t)))) ;

fig_1 = figure;
 % Extracting States [xL_dot, q_dot, reshape(R_dot, 9,1), vL_dot ,omega_dot, Omega_dot]   
subplot(2,2,1);
plot(t(ind),x(ind,1),'-b',t(ind),xLd(ind,1),':r',t(ind),x(ind,1)-x(ind,4),'-g');
grid on; title('x');legend('xL','x_dL','xQ');%axis equal;
xlabel('time');ylabel('x [m]');

subplot(2,2,2);
plot(t(ind),x(ind,2),'-b',t(ind),xLd(ind,2),':r',t(ind),x(ind,2)-x(ind,5),'-g');
grid on; title('y');legend('yL','y_dL','yQ');%axis equal;
xlabel('time');ylabel('y [m]');

subplot(2,2,3);
plot(t(ind),x(ind,3),'-b',t(ind),xLd(ind,3),':r',t(ind),x(ind,3)-x(ind,6),'-g');
grid on; title('z');legend('zL','z_dL','zQ');%axis equal;
xlabel('time');ylabel('z [m]');

subplot(2,2,4);
plot3(x(ind,1),x(ind,2),x(ind,3),'-g',xLd(ind,1),xLd(ind,2),xLd(ind,3),':r',x(ind,1)-x(ind,4),x(ind,2)-x(ind,5),x(ind,3)-x(ind,6),'-b');
grid on; title('trajectory');legend('traj_L','traj_d','traj_Q');%axis equal;
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
sgtitle({'Geometric Control and Differential Flatness of a Quadrotor UAV with a Cable-Suspended Load' , 'Koushil Sreenath, Taeyoung Lee, Vijay Kumar'});

if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
elseif ispc
    fig_1.WindowState = 'maximized';
    Image = getframe(fig_1);
    imwrite(Image.cdata, './figures/quad.jpg');
else
end

%% Animation
animate_3dquad_load(t, x, t(ind), xLd(ind,:));
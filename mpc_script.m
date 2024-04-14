% ex_twolink_drv.m

% Cleaning up variables from previous code
clearvars;
% Closing any open figures
close all;

% Defining model parameters

n = 3; %number of controlled joints
u0 = [0 0];
m = 3;
p = 10;
Qp = 300;
Qt = 500;
Ru=1;
t_max = 5;
dt_max = 0.01;

Q = Qp*eye(p);
Q(p,p) = Qt;
R = Ru*eye(m);

% Initial link angular position, in rad
theta_0 = zeros(n,1);
% Initial link velocity, in rad/s
theta_dot0 = zeros(n,1);
% Input voltage in V
%v = [500; 50];
%input
%u = [v(1)*km(1)/(r(1)*Ra(1));v(2)*km(2)/(r(2)*Ra(2))];
% Input torque frequency
%tau_freq = pi*[0.1; 0.2];
% Simulation stop time, in s
% Simulation solver maximum step time, in s


% Running the simulation
sim_out = sim("mpc.slx",'StopTime',num2str(t_max), ...
              'MaxStep',num2str(dt_max));

% Plotting the results
sim_out
% Time is the first column of the output matrix saved from the scope blocks
% in the model, the joint angles are the second and third columns
t_theta = sim_out.theta(:,1);
theta = sim_out.theta(:,2);
t_theta_dot = sim_out.theta_dot(:,1);
theta_dot = sim_out.theta_dot(:,2);

t_theta_d = sim_out.theta_d(:,1);
theta_d = sim_out.theta_d(:,2);
t_theta_d_dot = sim_out.theta_d_dot(:,1);
theta_d_dot = sim_out.theta_d_dot(:,2);

t_v = sim_out.V.Time;

v1 = sim_out.V.Data(1,1,:);

size(v1)
%t_tau_m = sim_out.tau_m(:,1);
%tau_m = sim_out.tau_m(:,2:3);

fig=figure(1);
subplot(3,1,1)
plot(t_theta, theta(:,1),LineWidth=1.5)
hold on
plot(t_theta, theta_d(:,1),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Pendulum')
grid on

subplot(3,1,2)
plot(t_theta, -theta(:,1)+theta_d(:,1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on

subplot(3,1,3)
plot(t_v, v1(1,:),LineWidth=1.5);
%ylim([0 max(v1(:,1))+100])
xlabel('Time [s]');
ylabel('Inputs')
grid on

exportgraphics(fig,"sliding_soft.png")

fig = figure(2);
plot(theta(:,1)-theta_d(:,1),theta_dot,LineWidth=1.5)
xlabel('\theta-\theta_d[s]')
ylabel('$\dot{\theta}$[rad]','Interpreter','latex')
title('Sliding Mode Surface')
grid on
exportgraphics(fig,"sliding_soft_surface.png")
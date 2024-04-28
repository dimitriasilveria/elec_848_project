%% SMC Contoller for UR5
% Cleaning up variables from previous code
clearvars;
% Closing any open figures
close all;
clc;

%% Defining model parameters

%Number of joints
n = 3;
%Gravity constant
g = 9.8;
%Sliding controller constants
mu_1 = 1;
mu_2 = 3;
m = 7;
p = 5;
k = 5;
l = 7;
%Controller constants
K1 = ones(n,1)*20;
K2 = ones(n,1)*2;
%Soft switching constants
delta = 0.5;
alpha = 4;
beta = 1;
h = 1/delta^(alpha/beta);
phi = 0.5;

% Simulation stop time, in s
t_max = 2;
% Simulation solver maximum step time, in s
dt_max = 0.01;
%Time vector
t = 0:dt_max:t_max;

%% Initialization
%Trajectory - Sinusoidal
q_ref = [1.75+1*sin(2*pi*1*t);-1.8+0.5*sin(2*pi*1*t);1.75+0.3*cos(2*pi*1*t)];
dq_ref = 2*pi*[1*cos(2*pi*1*t);0.5*cos(2*pi*1*t);-0.3*sin(2*pi*1*t)];
ddq_ref = 4*pi*pi*[-1*sin(2*pi*1*t);-0.5*sin(2*pi*1*t);-0.3*cos(2*pi*1*t)];

%Initialize system states
q = zeros(3,length(t));
% q(:,1) = [pi/2;-8*pi/9;5*pi/6]; %Initial offset
q(:,1) = q_ref(:,1);
dq = zeros(3,length(t));
ddq = zeros(3,length(t));

%Tracking Errors
e1 = zeros(3,length(t));
e2 = zeros(3,length(t));

%Control signals
U = zeros(3,length(t));
q_ddot_U = zeros(3,length(t));

%Cartesian Coordinates
p_Cart = zeros(3,length(t));
p_TrajCart = zeros(3,length(t));
e_Traj = zeros(1,length(t));
e_TrajXYZ = zeros(3,length(t));

%Torque disturbance 
% Tau_dist1 = 1*rand(n,length(t));
% Tau_dist = Tau_dist1 - mean(Tau_dist1);
Tau_dist = zeros(n,length(t));

%% SMC Controller
tic
for i = 1:length(t)   
    %Dynamics from MPC Model
    [T,Jacobi, M,C,G]  = fwdKIN(q(:,i), dq(:,i),g);
    %Control input - acceleration
    q_ddot_U(:,i)  = sliding_surface_to_acceleration(ddq_ref(:,i),e2(:,i),p, k, m, l, K1, K2, h, phi, e1(:,i), mu_2, mu_1);
    %Input torque
    U(:,i) = M*q_ddot_U(:,i) + C*dq(:,i) + G + Tau_dist(:,i); 
    %Saturation
    % U(:,i) = min(200, max(-200, U(:,i)));
    
    %Manipulator Dynamics
    ddq(:,i) = M\(U(:,i) - C*dq(:,i) - G);
    q(:,i+1) =  q(:,i) + dq(:,i)*dt_max + ddq(:,i)*dt_max*dt_max/2;
    dq(:,i+1) =  dq(:,i) + ddq(:,i)*dt_max;
    
    %Error updation
    e1(:,i+1) = q_ref(:,i) - q(:,i);  
    e2(:,i+1) = dq_ref(:,i) - dq(:,i);
    
    %Ideal trajectory
    T_traj  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);
    
    %Trajectory in cartesian space
    p_Cart(:,i) = real(T(1:3,4));
    p_TrajCart(:,i) = T_traj(1:3,4);
    e_Traj(:,i) = rms(p_TrajCart(:,i)-p_Cart(:,i));
    e_TrajXYZ(:,i) = (p_TrajCart(:,i)-p_Cart(:,i));
end
toc

%% PLOTS
%Joint 01 Trajectory
fig=figure(1);
subplot(3,1,1)
plot(t, q(1,1:end-1),LineWidth=1.5)
hold on
plot(t, q_ref(1,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 01')
grid on
subplot(3,1,2)
plot(t, e1(1,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
subplot(3,1,3)
plot(t, U(1,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
% exportgraphics(fig, "SMC_JointSpace_J1.png")

%Joint 02 Trajectory
fig=figure(2);
subplot(3,1,1)
plot(t, q(2,1:end-1),LineWidth=1.5)
hold on
plot(t, q_ref(2,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 02')
grid on
subplot(3,1,2)
plot(t, e1(2,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
subplot(3,1,3)
plot(t, U(2,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
% exportgraphics(fig, "SMC_JointSpace_J2.png")

%Joint 03 Trajectory
fig=figure(3);
subplot(3,1,1)
plot(t, q(3,1:end-1),LineWidth=1.5)
hold on
plot(t, q_ref(3,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 03')
grid on
subplot(3,1,2)
plot(t, e1(3,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
subplot(3,1,3)
plot(t, U(3,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
% exportgraphics(fig, "SMC_JointSpace_J3.png")

%Cartesian space trajectory
fig=figure(4);
plot3(p_TrajCart(1,2:end-1),p_TrajCart(2,2:end-1),p_TrajCart(3,2:end-1),'--',LineWidth=1.5);
hold on;grid on;
plot3(p_Cart(1,5:end-1),p_Cart(2,5:end-1),p_Cart(3,5:end-1),LineWidth=1.5);
hold off
xlabel('x[m]');
ylabel('y[m]');
zlabel('z[m]');
legend('Reference Trajectory','Tracking Trajectory','Location','best')
% exportgraphics(fig, "SMC_CartesianSpace.png")

%RMS Tracking error in cartesian space
fig=figure(5);
plot(t(:,1:end-1), e_Traj(:,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('RMS Tracking Error[m]')
grid on
title('Tracking Error')
% exportgraphics(fig, "SMC_TrackingError.png")

%Tracking error in cartesian space
fig = figure(6);
subplot(3,1,1)
plot(t(:,1:end-1), e_TrajXYZ(1,1:end-1),LineWidth=1.5)
xlabel('Time[s]')
ylabel('Error[m]')
title('X Axis Error')
grid on
subplot(3,1,2)
plot(t(:,1:end-1), e_TrajXYZ(2,1:end-1),LineWidth=1.5)
xlabel('Time[s]')
ylabel('Error[m]')
title('Y Axis Error')
grid on
subplot(3,1,3)
plot(t(:,1:end-1), e_TrajXYZ(3,1:end-1),LineWidth=1.5)
xlabel('Time[s]')
ylabel('Error[m]')
title('Z Axis Error')
grid on
% exportgraphics(fig, "SMC_TrackingErrorXYZ.png")
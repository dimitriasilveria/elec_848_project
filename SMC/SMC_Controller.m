%%
% Cleaning up variables from previous code
clearvars;
% Closing any open figures
close all;
clc

% Defining model parameters

n = 3; %number of controlled joints
mu_1 = 1;
mu_2 = 3;
m = 7;
p = 5;
% k = 5;
% l = 7;
k = 7;
l = 5;
K1 = ones(n,1)*20;
K2 = ones(n,1)*2;
delta = 0.5;
alpha = 4;
beta = 1;
h = 1/delta^(alpha/beta);
phi = 0.5;

% Simulation stop time, in s
t_max = 5;
% Simulation solver maximum step time, in s
dt_max = 0.01;
%Time vector
t = 0:dt_max:t_max;

g = 9.8;
%% Initialization

%Trajectory

q_ref = [2*sin(0.1*t);2*sin(0.5*t);2*sin(1*t)];
dq_ref = [2*0.1*cos(0.1*t);2*0.5*cos(0.5*t);2*1*cos(1*t)];
ddq_ref = [-2*0.1^2*sin(0.1*t);-2*0.5^2*sin(0.5*t);-2*1*sin(1*t)];
% qt = [pi/4;pi/6;-pi/4];
% q0 = [0;0;0];



%Initializa system states
q = zeros(3,length(t));
dq = zeros(3,length(t));
ddq = zeros(3,length(t));

%Tracking Errors
e1 = zeros(3,length(t));
e2 = zeros(3,length(t));

%Control signals
U = zeros(3,length(t));
q_ddot_U = zeros(3,length(t));
%% MPC Controller

for i = 1:length(t)
    
    %Dynamics from MPC Model
    [T,Jacobi, M,C,G]  = fwdKIN(q(:,i), dq(:,i),g);

    q_ddot_U(:,i)  = sliding_surface_to_acceleration(ddq_ref(:,i),e2(:,i),p, k, m, l, K1, K2, h, phi, e1(:,i), mu_2, mu_1);
    
    %Input torque
    U(:,i) = M*q_ddot_U(:,i) + C*dq(:,i) + G; %-tau_d;
    
    %Manipulator Dynamics
    ddq(:,i) = M\(U(:,i) - C*dq(:,i) - G);
    q(:,i+1) =  q(:,i) + dq(:,i)*dt_max + ddq(:,i)*dt_max*dt_max/2;
    dq(:,i+1) =  dq(:,i) + ddq(:,i)*dt_max;
    e1(:,i) = q_ref(:,i) - q(:,i);  
    % e2(:,i) = dq_ref(:,i) - dq(:,i); 

end


%% 

fig=figure(1);
subplot(3,3,1)
plot(t, q(1,1:end-1),LineWidth=1.5)
hold on
plot(t, q_ref(1,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 01')
grid on
% ylim([-2 +2])
subplot(3,3,2)
plot(t, q(2,1:end-1),LineWidth=1.5)
hold on
plot(t, q_ref(2,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 02')
grid on
% ylim([-2 +2])
subplot(3,3,3)
plot(t, q(3,1:end-1),LineWidth=1.5)
hold on
plot(t, q_ref(3,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 03')
grid on
% ylim([-2 +2])

subplot(3,3,4)
plot(t, e1(1,:),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,5)
plot(t, e1(2,:),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,6)
plot(t, e1(3,:),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])


subplot(3,3,7)
plot(t, U(1,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,8)
plot(t, U(2,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,9)
plot(t, U(3,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

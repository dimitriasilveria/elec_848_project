%% MPC Contoller for UR5

% Cleaning up variables from previous code
clearvars;
% Closing any open figures
close all;
%Clearing console
clc

%% Defining model parameters

%Gravitational constant [m/sec] 
g = 9.8;
%Number of controlled joints
n = 3;                                                                      
%Control horizon
m = 3;
%Prediction horizon
p = 10;
%Error weight matrix
% Qp = 300;
Qp = 300;
%Terminal weight matrix
% Qt = 500;
Qt = 500;
%Control weight matrix
% Ru=1;
Ru=100;

%Maximum Simulation time [sec]
t_max = 5;
%Sampling period [sec]
dt_max = 0.01;

%Time vector
t = 0:dt_max:t_max;

%Time
K0 = 100*ones(n,1);
K1 = 20*ones(n,1);

%% Initialization

%Trajectory
q_ref = [2*sin(0.1*t);2*sin(0.5*t);2*sin(1*t)];
dq_ref = [2*0.1*cos(0.1*t);2*0.5*cos(0.5*t);2*1*cos(1*t)];
ddq_ref = [-2*0.1^2*sin(0.1*t);-2*0.5^2*sin(0.5*t);-2*1*sin(1*t)];



%Initializa system states
q = zeros(3,length(t));
dq = zeros(3,length(t));
ddq = zeros(3,length(t));

%Tracking Errors
e1 = zeros(3,length(t));
e2 = zeros(3,length(t));

%Control signals
U_controller = zeros(3,length(t));
v = zeros(3,length(t));
%% MPC Controller

for i = 1:length(t)
    
    %Dynamics from MPC Model
    [T,Jacobi, M,C,G]  = fwdKIN(q(:,i), dq(:,i),g);
    % [T,Jacobi, M,C,G]  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);
           
    % Controller Input
    v(:,i) = ddq_ref(:,i)  +K0.*e1(:,i) + K1.*e2(:,i); 
    U_controller(:,i) = M*v(:,i) + C*dq(:,i) + G;
    
    %Manipulator Dynamics
    ddq(:,i) = M\(U_controller(:,i) - C*dq(:,i) - G);

    q(:,i+1) =  q(:,i) + dq(:,i)*dt_max + ddq(:,i)*dt_max*dt_max/2;
    dq(:,i+1) =  dq(:,i) + ddq(:,i)*dt_max;
    e1(:,i+1) = q_ref(:,i) - q(:,i);  
    e2(:,i+1) = dq_ref(:,i) - dq(:,i); 

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
plot(t, e1(1,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,5)
plot(t, e1(2,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,6)
plot(t, e1(3,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])


subplot(3,3,7)
plot(t, U_controller(1,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,8)
plot(t, U_controller(2,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,9)
plot(t, U_controller(3,:),LineWidth=1.5);
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

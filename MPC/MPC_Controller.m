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
Ru=0.1;

%Maximum Simulation time [sec]
t_max = 5;
%Sampling period [sec]
dt_max = 0.01;
%Predictive horizon [sec]
Np = p*dt_max;

%Time vector
t = 0:dt_max:t_max;
%Time vector for MPC
t_MPC = 0:dt_max:(t_max+Np);

%MPC Error Matrix
Q_p = Qp*ones(n*p-n,1);
Q_t = Qt*ones(n,1);
Q = [Q_p;Q_t];
Q = Q.*kron(eye(n),eye(p));

%MPC Control Matrix
R = Ru.*kron(eye(n),eye(m));


%% Initialization

%Trajectory
% q_ref = [2*sin(2*t_MPC);2*sin(2*t_MPC);2*sin(2*t_MPC)];
% dq_ref = [2*0.1*cos(0.1*t_MPC);2*0.5*cos(0.5*t_MPC);2*1*cos(1*t_MPC)];

% q_ref = [2*sin(0.1*t_MPC);2*sin(0.5*t_MPC);2*sin(1*t_MPC)];
% dq_ref = [2*0.1*cos(0.1*t_MPC);2*0.5*cos(0.5*t_MPC);2*1*cos(1*t_MPC)];

% q_ref = [t_MPC.*0.;t_MPC.*0;t_MPC.*0];
% %dq_ref = 0.1*ones(n,length(t_MPC));
% dq_ref = [0.1*zeros(1,length(t_MPC));zeros(2,length(t_MPC))];

% q_ref = [2*sin(2*t_MPC);2*sin(2*t_MPC);2*sin(2*t_MPC)];
% dq_ref = [2*0.1*cos(0.1*t_MPC);2*0.5*cos(0.5*t_MPC);2*1*cos(1*t_MPC)];

qt = [pi/4;pi/6;-pi/4];
q0 = [0;0;0];
traj_endPoint = [q0 qt];
t_traj_endPoint = [0 t_MPC(end)];
q_ref = zeros(3,length(t_MPC));

% for i = 1:n
%     q_ref(i,:) = interpn(t_traj_endPoint,traj_endPoint(i,:),t_MPC,'spline');
% end

for i = 1:n
    pol_traj = polyfit(t_traj_endPoint,traj_endPoint(i,:),5);
    q_ref(i,:) = polyval(pol_traj,t_MPC);
end

diff_q = diff(q_ref,1,2);
dq_ref = diff_q./dt_max;
dq_ref = [zeros(3,1) dq_ref];

%Initialize system states
q = zeros(3,length(t));
dq = zeros(3,length(t));
ddq = zeros(3,length(t));

%Tracking Errors
e1 = zeros(3,length(t));
e2 = zeros(3,length(t));

%Control signals
U_old = zeros(n*p,1);
U_controller = zeros(3,length(t));
tau = zeros(3,length(t));

%Cartesian Coordinates
p_Cart = zeros(3,length(t));
p_TrajCart = zeros(3,length(t));
e_Traj = zeros(1,length(t));
%% MPC Controller

for i = 2:length(t)
    
    %Dynamics from MPC Model
    [T,Jacobi, M,C,G]  = fwdKIN(q(:,i-1), dq(:,i-1),g);
    %[T,Jacobi, M,C,G]  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);

    Cqdot = C*dq(:,i-1);
    %Cqdot = C*dq_ref(:,i);

    [A,B,C,D] = state_space_matrices(M,Cqdot,G,n, dt_max);

    [S,W,V,L] = state_space_combine(p, m,n, Np,C,A,D,B);
    
    Y_ref = [];
    for j = 0:p-1
        % Y_ref = [Y_ref; q_ref(:,i+j);dq_ref(:,i+j)];
        Y_ref = [Y_ref; q_ref(:,i+j)];
    end
    
    x = [q(:,i-1);dq(:,i-1)];
    
    E = Y_ref - S*x - W*U_old - V;

    dU = MPC_optimization(E,W,R,L, Q);
    
    % dU_sat = min(8, max(-8, dU));
    dU_sat = dU;

    U = U_old + L*dU_sat;
    % U = min(20, max(-20, U));

    % U_old = U;
    U_old = repmat(U(1:3)',1,p)';
    
    % Controller Input
    U_controller(:,i) = U(1:3);
    
    % [T,Jacobi, M,C,G]  = fwdKIN(q(:,i-1), dq(:,i-1),g);
    % Cqdot = C*dq(:,i-1);
    tau(:,i) = M*U_controller(:,i) + Cqdot + G;
    
    %Mancipulator Dynamics
    ddq(:,i) = M\(tau(:,i) - Cqdot - G);

    q(:,i) =  q(:,i-1) + dq(:,i-1)*dt_max + ddq(:,i)*dt_max*dt_max/2;
    dq(:,i) =  dq(:,i-1) + ddq(:,i)*dt_max;
    e1(:,i) = q_ref(:,i) - q(:,i);  
    % e2(:,i) = dq_ref(:,i) - dq(:,i);
    
    T_traj  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);

    p_Cart(:,i-1) = T(1:3,4);
    p_TrajCart(:,i-1) = T_traj(1:3,4);
    e_Traj(:,i-1) = rms(p_TrajCart(:,i-1)-p_Cart(:,i-1));
end


%% 
fig=figure(1);
subplot(3,3,1)
plot(t, q(1,1:end),LineWidth=1.5)
hold on
plot(t_MPC, q_ref(1,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 01')
grid on
% ylim([-2 +2])
subplot(3,3,2)
plot(t, q(2,1:end),LineWidth=1.5)
hold on
plot(t_MPC, q_ref(2,:),'r--',LineWidth=1.5)
legend('Actual','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 02')
grid on
% ylim([-2 +2])
subplot(3,3,3)
plot(t, q(3,1:end),LineWidth=1.5)
hold on
plot(t_MPC, q_ref(3,:),'r--',LineWidth=1.5)
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
exportgraphics(fig, "MPC_JointSpace.png")

fig=figure(2);
plot3(p_TrajCart(1,2:end-1),p_TrajCart(2,2:end-1),p_TrajCart(3,2:end-1),'--',LineWidth=1.5);
hold on;grid on;
plot3(p_Cart(1,2:end-1),p_Cart(2,2:end-1),p_Cart(3,2:end-1),LineWidth=1.5);
hold off
xlabel('x[m]');
ylabel('y[m]');
zlabel('z[m]');
legend('Reference Trajecotry','Tracking Trajectory')
exportgraphics(fig, "MPC_CartesianSpace.png")

fig=figure(3);
plot(t(:,1:end-1), e_Traj(:,1:end-1),LineWidth=1.5);
xlabel('Time[s]')
ylabel('RMS Tracking Error[rad]')
grid on
exportgraphics(fig, "MPC_TrackingError.png")
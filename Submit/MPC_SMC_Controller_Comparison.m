%% MPC and SMC Contoller for UR5
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
%MPC Control horizon
m_MPC = 3;
%MPC Prediction horizon
p_MPC = 10;
%MPC Error weight matrix
Qp = 300;
%MPC Terminal weight matrix
Qt = 500;
%MPC Control weight matrix
Ru=1;
%Sliding controller constants
mu_1 = 1;
mu_2 = 3;
m = 7;
p = 5;
k = 5;
l = 7;
%SMC Controller constants
K1 = ones(n,1)*20;
K2 = ones(n,1)*2;
%SMC Soft switching constants
delta = 0.5;
alpha = 4;
beta = 1;
h = 1/delta^(alpha/beta);
phi = 0.5;

%Maximum Simulation time [sec]
t_max = 5;
%Sampling period [sec]
dt_max = 0.01;
%Predictive horizon [sec]
Np = p_MPC*dt_max;

%Time vector
t = 0:dt_max:t_max;
%Time vector for MPC
t_MPC = 0:dt_max:(t_max+Np);

%MPC Error Matrix
Q_p = Qp*ones(n*p_MPC-n,1);
Q_t = Qt*ones(n,1);
Q = [Q_p;Q_t];
Q = Q.*kron(eye(n),eye(p_MPC));

%MPC Control Matrix
R = Ru.*kron(eye(n),eye(m_MPC));

%Torque disturbance 
% Tau_dist1 = 1*rand(n,length(t_MPC));
% Tau_dist = Tau_dist1 - mean(Tau_dist1);
Tau_dist = zeros(n,length(t_MPC));

%% Initialization

%Trajectory 01 - 5th Order interpolation
qt = [pi/4;pi/6;-pi/4];
q0 = [0;0;0];
traj_endPoint = [q0 qt];
t_traj_endPoint = [0 t_MPC(end)];
q_ref = zeros(3,length(t_MPC));
for i = 1:n
    pol_traj = polyfit(t_traj_endPoint,traj_endPoint(i,:),5);
    q_ref(i,:) = polyval(pol_traj,t_MPC);
end
diff_q = diff(q_ref,1,2);
dq_ref = diff_q./dt_max;
dq_ref = [zeros(3,1) dq_ref];

diff_dq = diff(dq_ref,1,2);
ddq_ref = diff_dq./dt_max;
ddq_ref = [zeros(3,1) ddq_ref];

% %Trajectory 02 - Sinusoidal
% q_ref = [1.75+1*sin(2*pi*1*t_MPC);-1.8+0.5*sin(2*pi*1*t_MPC);1.75+0.3*cos(2*pi*1*t_MPC)];
% dq_ref = 2*pi*[1*cos(2*pi*1*t_MPC);0.5*cos(2*pi*1*t_MPC);-0.3*sin(2*pi*1*t_MPC)];
% ddq_ref = 4*pi*pi*[-1*sin(2*pi*1*t_MPC);-0.5*sin(2*pi*1*t_MPC);-0.3*cos(2*pi*1*t_MPC)];

%Initialize system states
q_MPC = zeros(3,length(t));
q_MPC(:,1) = q_ref(:,1);
% q_MPC(:,1) = [-pi/2;-pi/2;pi/2];
dq_MPC = zeros(3,length(t));
ddq_MPC = zeros(3,length(t));

q_SMC = zeros(3,length(t));
q_SMC(:,1) = q_ref(:,1);
% q_SMC(:,1) = [-pi/2;-pi/2;pi/2];
dq_SMC = zeros(3,length(t));
ddq_SMC = zeros(3,length(t));

%Tracking Errors
e1_MPC = zeros(3,length(t));
e1_MPC(:,1) = dq_ref(:,1) - q_MPC(:,1);
e2_MPC = zeros(3,length(t));

e1_SMC = zeros(3,length(t));
e1_SMC(:,1) = dq_ref(:,1) - q_MPC(:,1);
e2_SMC = zeros(3,length(t));

%Control signals
U_old_MPC = zeros(n*p_MPC,1);
U_controller_MPC = zeros(3,length(t));
tau_MPC = zeros(3,length(t));

U_SMC = zeros(3,length(t));
q_ddot_U_SMC = zeros(3,length(t));

%Cartesian Coordinates
p_Cart_MPC = zeros(3,length(t));
p_Cart_SMC = zeros(3,length(t));

p_TrajCart = zeros(3,length(t));

e_Traj_MPC = zeros(1,length(t));
e_Traj_SMC = zeros(1,length(t));

e_TrajXYZ_MPC = zeros(3,length(t));
e_TrajXYZ_SMC = zeros(3,length(t));

%% MPC Controller

for i = 2:length(t)
%% MPC
    %Dynamics from MPC Model
    [T,Jacobi, M,C,G]  = fwdKIN(q_MPC(:,i-1), dq_MPC(:,i-1),g);
    Cqdot = C*dq_MPC(:,i-1);
    %Discrete time state space matrices
    [A,B,C1,D] = state_space_matrices(M,Cqdot,G,n, dt_max);
    %Predictive Model
    [S,W,V,L] = state_space_combine(p_MPC, m_MPC,n, Np,C1,A,D,B);
    %Reference input over prediction horizon
    Y_ref = [];
    for j = 0:p_MPC-1
        Y_ref = [Y_ref; q_ref(:,i+j)];
    end
    %System state
    x = [q_MPC(:,i-1);dq_MPC(:,i-1)];
    %Error over prediction horizon
    E = Y_ref - S*x - W*U_old_MPC - V;
    
    %Incremental control input - QP solution
    dU = MPC_optimization(E,W,R,L, Q);
    
    %Saturation - Incremental
    % dU_sat = min(8, max(-8, dU));
    dU_sat = dU;
    %Control input - acceleration
    U = U_old_MPC + L*dU_sat;
    %Find matrix of old inputs
    U_old_MPC = repmat(U(1:3)',1,p_MPC)';
    % Controller Input
    U_controller_MPC(:,i) = U(1:3);
    
    %Torque input
    tau_MPC(:,i) = M*U_controller_MPC(:,i) + Cqdot + G + Tau_dist(:,i);
    %Toruqe saturation
    % tau_MPC(:,i) = min(100, max(-100, tau_MPC(:,i)));
    
    %Manipulator Dynamics
    ddq_MPC(:,i) = M\(tau_MPC(:,i) - Cqdot - G);
    q_MPC(:,i) =  q_MPC(:,i-1) + dq_MPC(:,i-1)*dt_max + ddq_MPC(:,i)*dt_max*dt_max/2;
    dq_MPC(:,i) =  dq_MPC(:,i-1) + ddq_MPC(:,i)*dt_max;
    e1_MPC(:,i) = q_ref(:,i) - q_MPC(:,i);  
    
    %Ideal trajectory
    T_traj  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);
    
    %Trajectory in cartesian space
    p_Cart_MPC(:,i-1) = T(1:3,4);
    p_TrajCart(:,i-1) = T_traj(1:3,4);
    e_Traj_MPC(:,i-1) = rms(p_TrajCart(:,i-1)-p_Cart_MPC(:,i-1));
    e_TrajXYZ_MPC(:,i-1) = (p_TrajCart(:,i-1)-p_Cart_MPC(:,i-1));

%% SMC
    %Dynamics from SMC Model   
    [T,Jacobi, M,C,G]  = fwdKIN(q_SMC(:,i-1), dq_SMC(:,i-1),g);
    %Control input - acceleration
    q_ddot_U_SMC(:,i-1)  = sliding_surface_to_acceleration(ddq_ref(:,i-1),e2_SMC(:,i-1),p, k, m, l, K1, K2, h, phi, e1_SMC(:,i-1), mu_2, mu_1);
    %Input torque
    U_SMC(:,i-1) = M*q_ddot_U_SMC(:,i-1) + C*dq_SMC(:,i-1) + G  + Tau_dist(:,i-1);
    %Saturation
    % U_SMC(:,i-1) = min(100, max(-100, real(U_SMC(:,i-1))));
   
    %Dynamics
    ddq_SMC(:,i-1) = M\(U_SMC(:,i-1) - C*dq_SMC(:,i-1) - G);
    dq_SMC(:,i) =  dq_SMC(:,i-1) + ddq_SMC(:,i-1)*dt_max;
    q_SMC(:,i) =  q_SMC(:,i-1) + dq_SMC(:,i)*dt_max + ddq_SMC(:,i-1)*dt_max*dt_max/2;
    
    %Error updation
    e1_SMC(:,i) = q_ref(:,i) - q_SMC(:,i);  
    e2_SMC(:,i) = dq_ref(:,i) - dq_SMC(:,i);
    
    %Trajectory in cartesian space
    p_Cart_SMC(:,i-1) = real(T(1:3,4));
    % p_TrajCart(:,i) = T_traj(1:3,4);
    e_Traj_SMC(:,i-1) = rms(p_TrajCart(:,i-1)-p_Cart_SMC(:,i-1));
    e_TrajXYZ_SMC(:,i-1) = (p_TrajCart(:,i-1)-p_Cart_SMC(:,i-1));

end

%% PLOTS
%Joint 01 Trajectory
fig=figure(1);
subplot(3,1,1)
plot(t, q_MPC(1,:),LineWidth=1.5)
hold on
plot(t, q_SMC(1,:),LineWidth=1.5)
plot(t_MPC, q_ref(1,:),'k--',LineWidth=1.5)
legend('MPC','SMC','Desired','Location','southeast')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 01')
grid on
subplot(3,1,2)
plot(t, e1_MPC(1,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(1,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
subplot(3,1,3)
plot(t, tau_MPC(1,:),LineWidth=1.5);
hold on
plot(t, U_SMC(1,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
% exportgraphics(fig, "COMB_JointSpace_J1.png")
% exportgraphics(fig, "COMB_JointSpace_J1_Offset.png")
% exportgraphics(fig, "COMB_JointSpace_J1_Tdist.png")

%Joint 02 Trajectory
fig=figure(2);
subplot(3,1,1)
plot(t, q_MPC(2,:),LineWidth=1.5)
hold on
plot(t, q_SMC(2,:),LineWidth=1.5)
plot(t_MPC, q_ref(2,:),'k--',LineWidth=1.5)
legend('MPC','SMC','Desired','Location','southeast')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 02')
grid on
subplot(3,1,2)
plot(t, e1_MPC(2,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(2,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
subplot(3,1,3)
plot(t, tau_MPC(2,:),LineWidth=1.5);
hold on
plot(t, U_SMC(2,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
% exportgraphics(fig, "COMB_JointSpace_J2.png")
% exportgraphics(fig, "COMB_JointSpace_J2_Offset.png")
% exportgraphics(fig, "COMB_JointSpace_J2_Tdist.png")

%Joint 03 Trajectory
fig=figure(3);
subplot(3,1,1)
plot(t, q_MPC(3,:),LineWidth=1.5)
hold on
plot(t, q_SMC(3,:),LineWidth=1.5)
plot(t_MPC, q_ref(3,:),'k--',LineWidth=1.5)
legend('MPC','SMC','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 03')
grid on
subplot(3,1,2)
plot(t, e1_MPC(3,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(3,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
subplot(3,1,3)
plot(t, tau_MPC(3,:),LineWidth=1.5);
hold on
plot(t, U_SMC(3,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
% exportgraphics(fig, "COMB_JointSpace_J3.png")
% exportgraphics(fig, "COMB_JointSpace_J3_Offset.png")
% exportgraphics(fig, "COMB_JointSpace_J3_Tdist.png")

%Cartesian space trajectory
fig=figure(4);
plot3(p_Cart_MPC(1,2:end-1),p_Cart_MPC(2,2:end-1),p_Cart_MPC(3,2:end-1),LineWidth=1.5);
hold on;grid on;
plot3(p_Cart_SMC(1,2:end-1),p_Cart_SMC(2,2:end-1),p_Cart_SMC(3,2:end-1),LineWidth=1.5);
plot3(p_TrajCart(1,2:end-1),p_TrajCart(2,2:end-1),p_TrajCart(3,2:end-1),'k--',LineWidth=1.5);
hold off
xlabel('x[m]');
ylabel('y[m]');
zlabel('z[m]');
legend('MPC Tracking Trajectory','SMC Tracking Trajectory','Reference Trajectory','Location','best')
% exportgraphics(fig, "COMP_CartesianSpace.png")
% exportgraphics(fig, "COMP_CartesianSpace_Tdist.png")
% exportgraphics(fig, "COMP_CartesianSpace_Offset.png")

%RMS Tracking error in cartesian space
fig=figure(5);
plot(t(:,1:end-1), e_Traj_MPC(:,1:end-1),LineWidth=1.5);
hold on
plot(t(:,1:end-1), e_Traj_SMC(:,1:end-1),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('RMS Tracking Error[m]')
grid on
title('Tracking Error')
% exportgraphics(fig, "COMP_TrackingError.png")
% exportgraphics(fig, "COMP_TrackingError_Tdist.png")
% exportgraphics(fig, "COMP_TrackingError_Offset.png")

%Tracking error in cartesian space
fig = figure(6);
subplot(3,1,1)
plot(t(:,1:end-1), e_TrajXYZ_MPC(1,1:end-1),LineWidth=1.5);
hold on
plot(t(:,1:end-1), e_TrajXYZ_SMC(1,1:end-1),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Error[m]')
title('X Axis Error')
grid on
subplot(3,1,2)
plot(t(:,1:end-1), e_TrajXYZ_MPC(2,1:end-1),LineWidth=1.5)
hold on
plot(t(:,1:end-1), e_TrajXYZ_SMC(2,1:end-1),LineWidth=1.5);
hold off
xlabel('Time[s]')
legend('MPC','SMC')
ylabel('Error[m]')
title('Y Axis Error')
grid on
subplot(3,1,3)
plot(t(:,1:end-1), e_TrajXYZ_MPC(3,1:end-1),LineWidth=1.5)
hold on
plot(t(:,1:end-1), e_TrajXYZ_SMC(3,1:end-1),LineWidth=1.5);
hold off
xlabel('Time[s]')
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Error[m]')
title('Z Axis Error')
grid on
% exportgraphics(fig, "COMP_TrackingErrorXYZ.png")
% exportgraphics(fig, "COMP_TrackingErrorXYZ_Tdist.png")
% exportgraphics(fig, "COMP_TrackingErrorXYZ_Offset.png")
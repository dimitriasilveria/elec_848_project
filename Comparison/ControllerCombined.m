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
m_MPC = 3;
%Prediction horizon
p_MPC = 10;
%Error weight matrix
% Qp = 300;
Qp = 300;
%Terminal weight matrix
% Qt = 500;
Qt = 500;
%Control weight matrix
% Ru=1;
Ru=0.1;

mu_1 = 1;
mu_2 = 3;
m = 7;
p = 5;
k = 5;
l = 7;
% k = 7;
% l = 5;
K1 = ones(n,1)*20;
K2 = ones(n,1)*2;
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

diff_dq = diff(dq_ref,1,2);
ddq_ref = diff_dq./dt_max;
ddq_ref = [zeros(3,1) ddq_ref];


%Initialize system states
q_MPC = zeros(3,length(t));
q_MPC(:,1) = q_ref(:,1);
dq_MPC = zeros(3,length(t));
ddq_MPC = zeros(3,length(t));

q_SMC = zeros(3,length(t));
q_SMC(:,1) = q_ref(:,1);
dq_SMC = zeros(3,length(t));
ddq_SMC = zeros(3,length(t));

%Tracking Errors
e1_MPC = zeros(3,length(t));
e2_MPC = zeros(3,length(t));

e1_SMC = zeros(3,length(t));
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


%% MPC Controller

for i = 2:length(t)
    
    %Dynamics from MPC Model
    [T,Jacobi, M,C,G]  = fwdKIN(q_MPC(:,i-1), dq_MPC(:,i-1),g);
    
    Cqdot = C*dq_MPC(:,i-1);
    
    [A,B,C1,D] = state_space_matrices(M,Cqdot,G,n, dt_max);

    [S,W,V,L] = state_space_combine(p_MPC, m_MPC,n, Np,C1,A,D,B);
    
    Y_ref = [];
    for j = 0:p_MPC-1
        Y_ref = [Y_ref; q_ref(:,i+j)];
    end
    
    x = [q_MPC(:,i-1);dq_MPC(:,i-1)];
    
    E = Y_ref - S*x - W*U_old_MPC - V;

    dU = MPC_optimization(E,W,R,L, Q);
    
    % dU_sat = min(8, max(-8, dU));
    dU_sat = dU;

    U = U_old_MPC + L*dU_sat;
    % U = min(20, max(-20, U));

    % U_old = U;
    U_old_MPC = repmat(U(1:3)',1,p_MPC)';
    
    % Controller Input
    U_controller_MPC(:,i) = U(1:3);
    tau_MPC(:,i) = M*U_controller_MPC(:,i) + Cqdot + G;
    
    
    %Manipulator Dynamics
    ddq_MPC(:,i) = M\(tau_MPC(:,i) - Cqdot - G);

    q_MPC(:,i) =  q_MPC(:,i-1) + dq_MPC(:,i-1)*dt_max + ddq_MPC(:,i)*dt_max*dt_max/2;
    dq_MPC(:,i) =  dq_MPC(:,i-1) + ddq_MPC(:,i)*dt_max;
    e1_MPC(:,i) = q_ref(:,i) - q_MPC(:,i);  
    % e2(:,i) = dq_ref(:,i) - dq(:,i);
    
    T_traj  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);

    p_Cart_MPC(:,i-1) = T(1:3,4);
    p_TrajCart(:,i-1) = T_traj(1:3,4);
    e_Traj_MPC(:,i-1) = rms(p_TrajCart(:,i-1)-p_Cart_MPC(:,i-1));
    
    %SMC   
    [T,Jacobi, M,C,G]  = fwdKIN(q_SMC(:,i-1), dq_SMC(:,i-1),g);

    q_ddot_U_SMC(:,i-1)  = sliding_surface_to_acceleration(ddq_ref(:,i-1),e2_SMC(:,i-1),p, k, m, l, K1, K2, h, phi, e1_SMC(:,i), mu_2, mu_1);
    
    %Input torque
    U_SMC(:,i-1) = M*q_ddot_U_SMC(:,i-1) + C*dq_SMC(:,i-1) + G; %-tau_d;
    %Dynamics
    ddq_SMC(:,i-1) = M\(U_SMC(:,i-1) - C*dq_SMC(:,i-1) - G);
    q_SMC(:,i) =  q_SMC(:,i-1) + dq_SMC(:,i-1)*dt_max + ddq_SMC(:,i-1)*dt_max*dt_max/2;
    dq_SMC(:,i) =  dq_SMC(:,i-1) + ddq_SMC(:,i-1)*dt_max;
    e1_SMC(:,i) = q_ref(:,i-1) - q_SMC(:,i-1);  
    e2_SMC(:,i) = dq_ref(:,i-1) - dq_SMC(:,i-1);

    p_Cart_SMC(:,i) = real(T(1:3,4));
    p_TrajCart(:,i) = T_traj(1:3,4);
    e_Traj_SMC(:,i) = rms(p_TrajCart(:,i)-p_Cart_SMC(:,i));


end


%% 
fig=figure(1);
subplot(3,3,1)
plot(t, q_MPC(1,:),LineWidth=1.5)
hold on
plot(t, q_SMC(1,:),LineWidth=1.5)
plot(t_MPC, q_ref(1,:),'k--',LineWidth=1.5)
legend('MPC','SMC','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 01')
grid on
% ylim([-2 +2])
subplot(3,3,2)
plot(t, q_MPC(2,:),LineWidth=1.5)
hold on
plot(t, q_SMC(2,:),LineWidth=1.5)
plot(t_MPC, q_ref(2,:),'k--',LineWidth=1.5)
legend('MPC','SMC','Desired')
hold off
xlabel('Time[s]')
ylabel('Angular Pos[rad]')
title('Joint 02')
grid on
% ylim([-2 +2])
subplot(3,3,3)
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
% ylim([-2 +2])

subplot(3,3,4)
plot(t, e1_MPC(1,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(1,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,5)
plot(t, e1_MPC(2,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(2,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,6)
plot(t, e1_MPC(3,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(3,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])


subplot(3,3,7)
plot(t, U_controller_MPC(1,:),LineWidth=1.5);
hold on
plot(t, U_SMC(1,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,8)
plot(t, U_controller_MPC(2,:),LineWidth=1.5);
hold on
plot(t, U_SMC(2,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,9)
plot(t, U_controller_MPC(3,:),LineWidth=1.5);
hold on
plot(t, U_SMC(3,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
exportgraphics(fig, "COMP_JointSpace.png")

fig=figure(2);
plot3(p_Cart_MPC(1,2:end-1),p_Cart_MPC(2,2:end-1),p_Cart_MPC(3,2:end-1),LineWidth=1.5);
hold on;grid on;
plot3(p_Cart_SMC(1,2:end-1),p_Cart_SMC(2,2:end-1),p_Cart_SMC(3,2:end-1),LineWidth=1.5);
plot3(p_TrajCart(1,2:end-1),p_TrajCart(2,2:end-1),p_TrajCart(3,2:end-1),'k--',LineWidth=1.5);
hold off
xlabel('x[m]');
ylabel('y[m]');
zlabel('z[m]');
legend('MPC Tracking Trajectory','SMC Tracking Trajectory','Reference Trajectory','Location','best')
exportgraphics(fig, "COMP_CartesianSpace.png")

fig=figure(3);
plot(t(:,1:end-1), e_Traj_MPC(:,1:end-1),LineWidth=1.5);
hold on
plot(t(:,1:end-1), e_Traj_SMC(:,1:end-1),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('RMS Tracking Error[rad]')
grid on
title('Tracking Error')
exportgraphics(fig, "COMP_TrackingError.png")
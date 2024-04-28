%% MPC Contoller for UR5

% Cleaning up variables from previous code
clearvars;
% Closing any open figures
close all;
%Clearing console
clc

%% Defining model parameters MPC

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
Ru=1;

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

%Trajectory 01
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

% %Trajectory 02
% q_ref = [1.75+1*sin(2*pi*1*t_MPC);-1.8+0.5*sin(2*pi*1*t_MPC);1.75+0.3*cos(2*pi*1*t_MPC)];
% dq_ref = 2*pi*[1*cos(2*pi*1*t_MPC);0.5*cos(2*pi*1*t_MPC);-0.3*sin(2*pi*1*t_MPC)];
% ddq_ref = 4*pi*pi*[-1*sin(2*pi*1*t_MPC);-0.5*sin(2*pi*1*t_MPC);-0.3*cos(2*pi*1*t_MPC)];

%Initialize system states
% q = zeros(3,length(t));
q(:,1) = [pi/2;-8*pi/9;5*pi/6];
q(:,1) = q_ref(:,1);
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
e_TrajXYZ = zeros(3,length(t));

%Torque disturbance 

% Tau_dist1 = 1*rand(n,length(t_MPC));
% Tau_dist = Tau_dist1 - mean(Tau_dist1);
Tau_dist = zeros(n,length(t_MPC));
% 

%% MPC Controller
tic
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
    tau(:,i) = M*U_controller(:,i) + Cqdot + G + Tau_dist(:,i);
    tau(:,i) = min(100, max(-100, tau(:,i)));

    %Manipulator Dynamics
    ddq(:,i) = M\(tau(:,i) - Cqdot - G);

    q(:,i) =  q(:,i-1) + dq(:,i-1)*dt_max + ddq(:,i)*dt_max*dt_max/2;
    dq(:,i) =  dq(:,i-1) + ddq(:,i)*dt_max;
    e1(:,i) = q_ref(:,i) - q(:,i);  
    % e2(:,i) = dq_ref(:,i) - dq(:,i);
    
    T_traj  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);

    p_Cart(:,i-1) = T(1:3,4);
    p_TrajCart(:,i-1) = T_traj(1:3,4);
    e_Traj(:,i-1) = rms(p_TrajCart(:,i-1)-p_Cart(:,i-1));
    e_TrajXYZ(:,i-1) = (p_TrajCart(:,i-1)-p_Cart(:,i-1));
end
toc

q_MPC = q;
e1_MPC = e1;
U_controller_MPC = tau;
e_Traj_MPC = e_Traj;
e_TrajXYZ_MPC = e_TrajXYZ;
p_Cart_MPC = p_Cart;

%% Defining model parameters MPC
n = 3; %number of controlled joints
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



% %
% mu_1 = 1;
% mu_2 = 3;
% m = 7;
% p = 5;
% k = 5;
% l = 7;
% % k = 7;
% % l = 5;
% K1 = ones(n,1)*20;
% K2 = ones(n,1)*2;
% delta = 0.5;
% alpha = 4;
% beta = 1;
% h = 1/delta^(alpha/beta);
% phi = 0.5;


%Initialize system states
q = zeros(3,length(t));
q(:,1) = [pi/2;-8*pi/9;5*pi/6];
% q(:,1) = q_ref(:,1);
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

%% SMC Controller
tic
for i = 1:length(t)
    
    %Dynamics from MPC Model
    [T,Jacobi, M,C,G]  = fwdKIN(q(:,i), dq(:,i),g);

    q_ddot_U(:,i)  = sliding_surface_to_acceleration(ddq_ref(:,i),e2(:,i),p, k, m, l, K1, K2, h, phi, e1(:,i), mu_2, mu_1);
    % q_ddot_U(:,i) = abs(q_ddot_U(:,i));
    %Input torque
    U(:,i) = M*q_ddot_U(:,i) + C*dq(:,i) + G + Tau_dist(:,i); %-tau_d;
    
    U(:,i) = min(100, max(-100, real(U(:,i))));

    % disp(M);
    % disp(C);
    % disp(G);
    % disp(Jacobi);
    %Manipulator Dynamics
    ddq(:,i) = M\(U(:,i) - C*dq(:,i) - G);
    q(:,i+1) =  q(:,i) + dq(:,i)*dt_max + ddq(:,i)*dt_max*dt_max/2;
    dq(:,i+1) =  dq(:,i) + ddq(:,i)*dt_max;
    e1(:,i+1) = q_ref(:,i) - q(:,i);  
    e2(:,i+1) = dq_ref(:,i) - dq(:,i);

    T_traj  = fwdKIN(q_ref(:,i), dq_ref(:,i),g);

    p_Cart(:,i) = real(T(1:3,4));
    p_TrajCart(:,i) = T_traj(1:3,4);
    e_Traj(:,i) = rms(p_TrajCart(:,i)-p_Cart(:,i));
    e_TrajXYZ(:,i) = (p_TrajCart(:,i)-p_Cart(:,i));
end
toc

q_SMC = q;
e1_SMC = e1;
U_controller_SMC = U;
e_Traj_SMC = e_Traj;
e_TrajXYZ_SMC = e_TrajXYZ;
p_Cart_SMC = p_Cart;


%% 
fig=figure(1);
subplot(3,3,1)
plot(t, q_MPC(1,:),LineWidth=1.5)
hold on
plot(t, q_SMC(1,1:end-1),LineWidth=1.5)
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
plot(t, q_SMC(2,1:end-1),LineWidth=1.5)
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
plot(t, q_SMC(3,1:end-1),LineWidth=1.5)
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
plot(t, e1_SMC(1,1:end-1),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,5)
plot(t, e1_MPC(2,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(2,1:end-1),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])
subplot(3,3,6)
plot(t, e1_MPC(3,:),LineWidth=1.5);
hold on
plot(t, e1_SMC(3,1:end-1),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time[s]')
ylabel('Angular Error[rad]')
grid on
% ylim([-2 +2])


subplot(3,3,7)
plot(t, tau(1,:),LineWidth=1.5);
hold on
plot(t, U_controller_SMC(1,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,8)
plot(t, tau(2,:),LineWidth=1.5);
hold on
plot(t, U_controller_SMC(2,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on

subplot(3,3,9)
plot(t, tau(3,:),LineWidth=1.5);
hold on
plot(t, U_controller_SMC(3,:),LineWidth=1.5);
hold off
legend('MPC','SMC')
xlabel('Time [s]');
ylabel('Torque [N.m]')
grid on
% exportgraphics(fig, "COMP_JointSpace.png")
% exportgraphics(fig, "COMP_JointSpace_Tdist.png")
% exportgraphics(fig, "COMP_JointSpace_Offset.png")

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
% exportgraphics(fig, "COMP_CartesianSpace.png")
% exportgraphics(fig, "COMP_CartesianSpace_Tdist.png")
% exportgraphics(fig, "COMP_CartesianSpace_Offset.png")

fig=figure(3);
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

fig = figure(4);
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
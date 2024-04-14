% Cleaning up variables from previous code
clearvars;
% Closing any open figures
close all;
clc

% Defining model parameters
g = 9.8;
n = 3; %number of controlled joints

m = 3;
p = 10;
Qp = 300;
Qt = 500;
Ru=1;
t_max = 5;
dt_max = 0.01;
Np = p*dt_max;

t = 0:dt_max:t_max;
t_MPC = 0:dt_max:(t_max+Np);

Q_p = Qp*ones(n*p-n,1);
Q_t = Qt*ones(n,1);
Q = [Q_p;Q_t];
Q = Q.*kron(eye(n),eye(p));

R = Ru.*kron(eye(n),eye(m));

% Initial link angular position, in rad
theta_0 = zeros(n,1);
% Initial link velocity, in rad/s
theta_dot0 = zeros(n,1);
% Input voltage in V

%% Dummy Feedback

q_ref = [sin(0.1*t_MPC);sin(0.5*t_MPC);sin(1*t_MPC)];
dq_ref = [0.1*cos(0.1*t_MPC);0.5*cos(0.5*t_MPC);1*cos(1*t_MPC)];

q = zeros(3,length(t));
dq = zeros(3,length(t));

U_old = zeros(n*p,1);

%% MPC Controller

for i = 1:length(t)
    
    %Dynamics from MPC Model
    
    [T,Jacobi, M,C,G]    = fwdKIN(q(:,i), dq(:,i),g);
    Cqdot = C*dq_ref(:,i);
    [A,B,C,D] = state_space_matrices(M,Cqdot,G,n, dt_max);

    [S,W,V,L] = state_space_combine(p, m,n, Np,C,A,D,B);
    
    Y_ref = [];
    for j = 0:p-1
        % Y_ref = [Y_ref; q_ref(:,i+j);dq_ref(:,i+j)];
        Y_ref = [Y_ref; q_ref(:,i+j)];
    end
    
    x = [q(:,i);dq(:,i)];
    
    E = Y_ref - S*x - W*U_old - V;

    dU = MPC_optimization(E,W,R,L, Q);

    U = U_old + L*dU;

    U_old = U;
    
    % Controller Input
    U_controller = U(1:3);

    %Manipulator Dynamics
    ddq = M\(U_controller - Cqdot - G);
    
    dq = 


end

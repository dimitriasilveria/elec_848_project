function [T,Jacobi, M,C,G]    = fwdKIN_6Joint(q, dq, g)

%Position and Velocity vectors
q_1 = q(1);
q_2 = q(2);
q_3 = q(3);
q_4 = q(4);
q_5 = q(5);
q_6 = q(6);
dq_1 = dq(1);
dq_2 = dq(2);
dq_3 = dq(3);

% TRANSLATION MATRICES AND FORWARD KINEMATICS 
% (Position of EE from joint space)
T = 

% JACOBIAN MATRIX OF END EFFECTOR
Jacobi = 
%Robots Inertia (Mass Matrix)
M =                                                
% CORIOLIS and CENTRIFUGAL MATRIX
C = 
%Gravity
G =
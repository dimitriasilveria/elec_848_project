function [A,B,C,D] = state_space_matrices(M,Cqdot,G,n, dt_max)
%This function implements the discrete time state space equation for the
%system

I_n = eye(n);
O_n = zeros(n);
k = inv(M);
g = -inv(M)*(Cqdot + G);

A = [I_n dt_max*I_n;O_n I_n];
B = [(dt_max^2/2)*k;dt_max*k];
C = [I_n O_n];
D = [(dt_max^2/2)*g;dt_max*g];
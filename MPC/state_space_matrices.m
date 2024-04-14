function [A,B,C,D] = state_space_matrices(M,Cqdot,G,n, dt_max)
%n is the number of joints 
I_n = eye(n);
O_n = zeros(n);
k = inv(M);
A = [I_n dt_max*I_n;O_n I_n];
B = [(dt_max^2/2)*k;dt_max*k];
C = [I_n O_n];
g = -inv(M)*(Cqdot + G); %here I am also assuming C = C(q,q_dot)*q_dot

D = [(dt_max^2/2)*g;dt_max*g];
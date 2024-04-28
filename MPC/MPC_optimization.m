function dU = MPC_optimization(E,W,R,L, Q)
% QP optimization of model parameters to find incremental control input

H = 2*(L'*W'*Q'*Q*W*L + R'*R);
H=(H+H')/2;
f = -2*L'*W'*Q'*Q*E;

dU = quadprog(H,f);


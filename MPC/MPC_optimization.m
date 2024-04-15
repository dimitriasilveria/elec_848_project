function dU = MPC_optimization(E,W,R,L, Q)

% H = 0.5*(L'*W'*Q'*Q*W*L + R'*R);
H = 2*(L'*W'*Q'*Q*W*L + R'*R);
H=(H+H')/2;
f = -2*L'*W'*Q'*Q*E;

dU = quadprog(H,f);
% dU = quadprog(H,f,[],[],[],[],-8,8);

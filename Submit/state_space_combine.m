function [S,W,V,L] = state_space_combine(p, m,n, Np,C,A,D,B)
%This function finds the system states over the prediction horizon and
%control input over the control horizon

S = [];
W = [];
V = [];
L = [];
C1 = [];
C2 = [];

for i = 1:p
    S = [S;C*A^i];
end
l_v = C*D.*0;
k = p;
for i = 1:p
   if i == 1
       C1 = [eye(n)];
       C2 = [zeros(n)];
   else
       C1 = [C1; eye(n)];
       C2 = [C2; eye(n)];
   end

   l_v = l_v + C*A^(p-k)*D;
   V = [V;l_v];   
   k = k - 1;
   l = [];
   for j = 1:p
       if j<=i
           l = [l C*A^(i-j)*B];
       else
           l = [l 0.*C*A*B];
       end
  
   end
   W = [W;l];
end
for i=1:m
    if i == m
        L = [L C1];
    else
        L = [L C2];
    end

end
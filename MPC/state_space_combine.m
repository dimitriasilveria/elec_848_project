function [S,W,V,L] = state_space_combine(p, m,n, Np,C,A,D,B)

S = [];
W = [];
V = [];
L = [];
L1 = [];
L2 = [];

for i = 1:p
    S = [S;C*A.^i];
end
l_v = C*D*0;
k = p;
for i = 1:p
   if i < p
       L1 = [L1; eye(n)];
   else
       L1 = [zeros(n); L1];
   end
   L2 = [L2; eye(n)];
   l_v = l_v + C*A.^(p-k)*D;
   V = [V;l_v];   
   k = k - 1;
   l = [];
   for j = 1:p
       if j<=i
           l = [l C*A.^(i-j)*B];
       else
           l = [l 0*C*A*B];
       end
  
   end
   W = [W;l];
end
for i=1:m
    if i == m
        L = [L L1];
    else
        L = [L L2];
    end

end



%Np = P*dt_max;
%Nm = m*d_t_max;
% S = [];
% W = [];
% V = [];
% L = [];
% 
% for i = 1:p
%     S = [S;C*A.^i];
% end
% 
% l_v = C*D.*0;
% k = p;
% 
% for i = 1:p
%    l_v = l_v + C*A.^(p-k)*D;
%    V = [V;l_v];   
%    k = k - 1;
%    l = [];
%    for j = 1:p
%        if j<=i
%            l = [l C*A.^(i-j)*B];
%        else
%            l = [l 0.*C*A*B];
%        end
% 
%    end
%    W = [W;l];
% end


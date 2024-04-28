function q_ddot  = sliding_surface_to_acceleration(q_d_ddot,e_dot,p, k, m, l, K1, K2, h, phi, e, mu_2, mu_1)
%This function provides the acceleration control input for the SMC
%controller

%Sliding surface
s =e_dot + mu_1*e.^(m/p)+mu_2*e.^(k/l);
sat = e.*0;

%Soft switching
for i =1:length(sat)

    if abs(s(i))>=phi
        sat(i) = sign(s(i));
    else
        sat(i)=h*s(i);
    end
end

%Reaching law
s_dot = -K1.*s-K2.*sat;

%Control input
for i =1:length(sat)
    if e(i)== 0
        q_ddot(i) = q_d_ddot(i);
    else
        q_ddot(i) = q_d_ddot(i) + mu_1*(m/p)*e(i).^((m/p)-1).*e_dot(i) + mu_2*(k/l)*e(i).^((k/l)-1).*e_dot(i) - s_dot(i);
    end
end

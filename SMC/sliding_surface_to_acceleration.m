function q_ddot  = sliding_surface_to_acceleration(q_d_ddot,e_dot,p, k, m, l, K1, K2, h, phi, e, mu_2, mu_1)
s =e_dot + mu_1*e.^(m/p)+mu_2*e.^(k/l);
sat = e.*0;

for i =1:length(sat)

    if abs(s(i))>=phi
        sat(i) = sign(s(i));
    else
        sat(i)=h*s(i);
    end
end
s_dot = -K1.*s-K2.*sat;

q_ddot = q_d_ddot + mu_1*(m/p)*e.^((m/p)-1).*e_dot + mu_2*(k/l)*e.^((k/l)-1).*e_dot - s_dot;
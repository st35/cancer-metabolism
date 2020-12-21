function r_G6PD = Glucose_6_Phosphate_Dehydrogenase(C_G6P, C_NADP, C_6PGDL, C_NADPH, Rate)
    k1 = 1.1e8; k3 = 0.26e8; k8 = 11.0e8; k10 = 14.0e8;
    k2 = 0.87e3; k4 = 0.30e3; k5 = 0.75e3; k6 = 2.0e3; k7 = 220.0e3; k9 = 10.0e3;
    Et = 1.0e-3;

    Et = Et*2.0;

    Et = Et*Rate;

    k1 = k1*3600.0*1e-3; k3 = k3*3600.0*1e-3; k8 = k8*3600.0*1e-3; k10 = k10*3600.0*1e-3;
    k2 = k2*3600.0; k4 = k4*3600.0; k5 = k5*3600.0; k6 = k6*3600.0; k7 = k7*3600.0; k9 = k9*3600.0;

    N1 = k1*k3*k5*k7*k9;
    N2 = k2*k4*k6*k8*k10;
    D1 = (k2*k9)*(k4*k6 + k5*k6 + k5*k7);
    D2 = (k1*k9)*(k4*k6 + k5*k6 + k5*k7);
    D3 = k3*k5*k7*k9;
    D4 = k2*k4*k6*k8;
    D5 = (k2*k10)*(k4*k6 + k5*k6 + k5*k7);
    D6 = (k1*k3)*(k5*k7 + k5*k9 + k6*k9 + k7*k9);
    D7 = k1*k4*k6*k8;
    D8 = k3*k5*k7*k10;
    D9 = (k8*k10)*(k2*k4 + k2*k5 + k2*k6 + k4*k6);
    D10 = (k1*k3*k8)*(k5 + k6);
    D11 = (k3*k8*k10)*(k5 + k6);

    A = C_NADP;
    B = C_G6P;
    P = C_6PGDL;
    Q = C_NADPH;

    t0 = Et*(N1*A*B - N2*P*Q);
    t1 = D1 + D2*A + D3*B + D4*P + D5*Q + D6*A*B + D7*A*P + D8*B*Q + D9*P*Q + D10*A*B*P + D11*B*P*Q;

    r_G6PD = t0 / t1;
end
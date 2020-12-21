function r_PDH = Pyruvate_Dehydrogenase(M_PYR, M_COASH, M_NAD, M_CO2, M_ACCOA, M_NADH, Rate)
    V_PDH_mf = (1.22e-1)*1e3*3600.0; V_PDH_mr = 0.0;
    K_PDH_A = 38.3e-3; K_PDH_B = 9.9e-3; K_PDH_C = 60.7e-3;
    K_PDH_ACCOA = 40.2e-3; K_PDH_NADH = 40.0e-3;

    V_PDH_mf = V_PDH_mf*0.5e2;
    V_PDH_mr = V_PDH_mr*0.5e2;

    V_PDH_mf = V_PDH_mf*Rate;
    V_PDH_mr = V_PDH_mr*Rate;

    A = M_PYR; B = M_COASH; C = M_NAD;
    P = M_CO2; Q = M_ACCOA; R = M_NADH;

    alpha_i1 = 1.0 + (M_ACCOA / K_PDH_ACCOA);
    alpha_i2 = (M_NADH / K_PDH_NADH);

    t0 = (V_PDH_mf*A*B*C) - (V_PDH_mr*P*Q*R);
    t1 = K_PDH_C*A*B + K_PDH_B*alpha_i1*A*C + K_PDH_A*B*C + A*B*C + alpha_i2;

    r_PDH = t0 / t1;
end
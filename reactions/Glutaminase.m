function r_GLS = Glutaminase(M_Q, M_GLU)
    V_GLS_mf = 39.20*60.0;
    K_GLS_Q = 6.0;

    K_GLS_GLU = 5.0e-2;

    t0 = (V_GLS_mf*M_Q / K_GLS_Q);
    t1 = 1.0 + (M_Q / K_GLS_Q) + (M_GLU / K_GLS_GLU);

    r_GLS = t0 / t1;
end
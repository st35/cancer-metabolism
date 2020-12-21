function r_CIT = Citrate_Synthase(M_ACCOA, M_OAA, M_CIT, M_COASH, M_ATP, M_ADP, M_AMP, M_SCOA, Rate)
    V_CIT_mf = 11.6*1e3*3600.0; V_CIT_mr = V_CIT_mf*1e2;
    K_CIT_A = 4.0e-3; K_CIT_B = 14.0e-3;
    K_CIT_ia = 3.33e-3;
    K_CIT_CIT = 1600.0e-3; K_CIT_ATP = 900.0e-3; K_CIT_ADP = 1800.0e-3; K_CIT_AMP = 6000.0e-3; K_CIT_COASH = 67.0e-3; K_CIT_SCOA = 140.0e-3;

    V_CIT_mf = V_CIT_mf*Rate;
    V_CIT_mr = V_CIT_mr*Rate;

    A = M_OAA; B = M_ACCOA; P = M_COASH; Q = M_CIT;

    alpha_i1 = 1.0 + (M_CIT / K_CIT_CIT);
    alpha_i2 = 1.0 + (M_ATP / K_CIT_ATP) + (M_ADP / K_CIT_ADP) + (M_AMP / K_CIT_AMP) + (M_COASH / K_CIT_COASH) + (M_SCOA / K_CIT_SCOA);

    t0 = (V_CIT_mf*A*B) - (V_CIT_mr*P*Q);
    t1 = 1.0 + K_CIT_ia*K_CIT_B*alpha_i1 + K_CIT_A*alpha_i1*B + K_CIT_B*alpha_i2*A + A*B;

    r_CIT = t0 / t1;
end
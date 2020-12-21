function r_SDH = Succinate_Dehydrogenase(M_SUC, M_COQ, M_FUM, M_COQH2, M_OAA, Rate)
    V_SDH_mf = 1.067*1e4*60.0; V_SDH_mr = 1.174*1e3*60.0;
    K_SDH_SUC = 467.0e-3; K_SDH_COQ = 480.0e-3; K_SDH_FUM = 1200.0e-3; K_SDH_COQH2 = 2.45e-3;
    K_SDH_SUCa = 450.0e-3; K_SDH_FUMa = 375.0e-3;
    K_SDH_OAAi = 1.5e-3; K_SDH_SUCi = 120.0e-3; K_SDH_FUMi = 1275.0e-3;

    V_SDH_mf = V_SDH_mf*Rate;
    V_SDH_mr = V_SDH_mr*Rate;

    t0 = ((V_SDH_mf*M_SUC*M_COQ) / (K_SDH_SUC*K_SDH_COQ));
    t1 = ((V_SDH_mr*M_FUM*M_COQH2) / (K_SDH_FUM*K_SDH_COQH2));

    t2 = 1.0 + (M_SUC / K_SDH_SUC) + (M_COQ / K_SDH_COQ);
    t3 = ((M_SUC*M_COQ) / (K_SDH_SUC*K_SDH_COQ));
    t4 = (M_FUM / K_SDH_FUM) + (M_COQH2 / K_SDH_COQH2);
    t5 = ((M_FUM*M_COQH2) / (K_SDH_FUM*K_SDH_COQH2));
    
    alpha_i1 = (1.0 + (M_OAA / K_SDH_OAAi) + (M_SUC / K_SDH_SUCi) + (M_FUM / K_SDH_FUMi));
    alpha_i2 = (1.0 + (M_SUC / K_SDH_SUCa) + (M_FUM / K_SDH_FUMa));
    alpha_i = alpha_i1 / alpha_i2;

    r_SDH = (t0 - t1) / (t2 + t3 + t4 + t5 + alpha_i);
end
function r_ACLY = ATP_Citrate_Lyase(C_CIT, C_MgATP, C_COA, C_ACCOA, C_MgADP, C_OAA, Rate)
    V_ACLY_mf = 1.8e4;
    K_ACLY_CIT = 7.0e-3; K_ACLY_MgATP = 150.0e-3; K_ACLY_COA = 12.0e-3;
    V_ACLY_mr = V_ACLY_mf / 1e2;
    K_ACLY_ACCOA = 5.0e-3; K_ACLY_MgADP = 0.46e-2; K_ACLY_OAA = 1e-3;

    V_ACLY_mf = V_ACLY_mf*Rate;
    V_ACLY_mr = V_ACLY_mr*Rate;

    t0 = ((V_ACLY_mf*C_CIT*C_MgATP*C_COA) / (K_ACLY_CIT*K_ACLY_MgATP*K_ACLY_COA));
    t1 = 1.0 + (C_CIT / K_ACLY_CIT) + (C_MgATP / K_ACLY_MgATP) + (C_COA / K_ACLY_COA);
    t2 = ((C_CIT*C_MgATP) / (K_ACLY_CIT*K_ACLY_MgATP));
    t3 = ((C_CIT*C_MgATP*C_COA) / (K_ACLY_CIT*K_ACLY_MgATP*K_ACLY_COA));
    t4 = ((V_ACLY_mr*C_ACCOA*C_MgADP*C_OAA) / (K_ACLY_ACCOA*K_ACLY_MgADP*K_ACLY_OAA));
    t5 = (C_ACCOA / K_ACLY_ACCOA) + (C_MgADP / K_ACLY_MgADP) + (C_OAA / K_ACLY_OAA);

    r_ACLY = (t0 - t4) / (t1 + t2 + t3 + t5);
end
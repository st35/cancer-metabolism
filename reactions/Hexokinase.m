function r_HK = Hexokinase(C_GLC, C_MgATP, C_G6P, C_MgADP, C_G16BP, C_23BPG, C_GSH, Rate)
    V_HK_mf = 6.38e3; V_HK_mr = 41.0;
    K_HK_GLC = 0.1; K_HK_MgATP = 1.0; K_HK_G6P = 0.047; K_HK_MgADP = 1.0;
    K_HK_G6P_i = 0.047; K_HK_G16BP = 0.03; K_HK_23BPG = 4.0; K_HK_GSH = 3.0;

    V_HK_mf = V_HK_mf*Rate;
    V_HK_mr = V_HK_mr*Rate;

    t0 = ((V_HK_mf*C_GLC*C_MgATP) / (K_HK_GLC*K_HK_MgATP)) - ((V_HK_mr*C_G6P*C_MgADP) / (K_HK_G6P*K_HK_MgADP));
    t1 = 1.0 + (C_GLC / K_HK_GLC) + (C_MgATP / K_HK_MgATP) + ((C_GLC*C_MgATP) / (K_HK_GLC*K_HK_MgATP));
    t2 = (C_G6P / K_HK_G6P) + (C_MgADP / K_HK_MgADP) + ((C_G6P*C_MgADP) / (K_HK_G6P*K_HK_MgADP));
    i0 = (C_GLC*C_G6P) / (K_HK_GLC*K_HK_G6P_i);
    i1 = (C_GLC*C_G16BP) / (K_HK_GLC*K_HK_G16BP);
    i2 = (C_GLC*C_23BPG) / (K_HK_GLC*K_HK_23BPG);
    i3 = (C_GLC*C_GSH) / (K_HK_GLC*K_HK_GSH);

    r_HK = (t0) / (t1 + t2 + i0 + i1 + i2 + i3);
end
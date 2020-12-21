function r_GAPDH = Glyceraldehyde_3_Phosphate_Dehydrogenase(C_GAP, C_NAD, C_Pi, C_13BPG, C_NADH, C_H)
    V_GAPDH_mf = 5.317e3; V_GAPDH_mr = 3.919e3;
    K_GAPDH_GAP = 0.095; K_GAPDH_NAD = 0.045; K_GAPDH_Pi = 3.16;
    K_GAPDH_13BPG = 0.0067; K_GAPDH_NADH = 0.0033;
    K_GAPDH_Pi_i = 3.16; K_GAPDH_GAP_i = 1.59e-16; K_GAPDH_13BPG_i = 1.52e-18; K_GAPDH_NADH_i = 0.01; K_GAPDH_NAD_i = 0.045;
    K_GAPDH_GAP_dash_i = 0.031; K_GAPDH_13BPG_dash_i = 0.001;

    V_GAPDH_mf = V_GAPDH_mf*1e2;
    V_GAPDH_mr = V_GAPDH_mr*1e2;

    t0 = ((V_GAPDH_mf*C_GAP*C_NAD*C_Pi) / (K_GAPDH_GAP_i*K_GAPDH_NAD*K_GAPDH_Pi_i));
    t1 = ((V_GAPDH_mr*C_13BPG*C_NADH*C_H) / (K_GAPDH_13BPG_i*K_GAPDH_NADH));
    t2 = (C_GAP / K_GAPDH_GAP_i)*(1.0 + (C_GAP / K_GAPDH_GAP_dash_i));
    t3 = (C_13BPG / K_GAPDH_13BPG_i)*(1.0 + (C_GAP / K_GAPDH_GAP_dash_i));
    t4 = (K_GAPDH_13BPG*C_NADH*C_H) / (K_GAPDH_13BPG_i*K_GAPDH_NADH);
    t5 = (K_GAPDH_GAP*C_NAD*C_Pi) / (K_GAPDH_NAD*K_GAPDH_Pi_i*K_GAPDH_GAP_i);
    t6 = (C_NAD*C_GAP) / (K_GAPDH_NAD_i*K_GAPDH_GAP_i);
    t7 = ((C_Pi*C_GAP) / (K_GAPDH_Pi_i*K_GAPDH_GAP_i))*(1.0 + (C_GAP / K_GAPDH_GAP_dash_i));
    t8 = (C_NAD*C_13BPG) / (K_GAPDH_NAD_i*K_GAPDH_13BPG_i);
    t9 = (K_GAPDH_13BPG*C_Pi*C_NADH*C_H) / (K_GAPDH_Pi_i*K_GAPDH_13BPG_i*K_GAPDH_NADH);
    t10 = (C_GAP*C_NADH*C_H) / (K_GAPDH_GAP_i*K_GAPDH_NADH_i);
    t11 = (C_13BPG*C_NADH*C_H) / (K_GAPDH_13BPG_i*K_GAPDH_NADH);
    t12 = (C_NAD*C_Pi*C_GAP) / (K_GAPDH_NAD*K_GAPDH_Pi_i*K_GAPDH_GAP_i);
    t13 = (K_GAPDH_GAP*C_NAD*C_Pi*C_13BPG) / (K_GAPDH_GAP_i*K_GAPDH_NAD*K_GAPDH_Pi_i*K_GAPDH_13BPG_dash_i);
    t14 = (C_Pi*C_GAP*C_NADH*C_H) / (K_GAPDH_Pi_i*K_GAPDH_GAP_i*K_GAPDH_NADH_i);
    t15 = (K_GAPDH_13BPG*C_Pi*C_13BPG*C_NADH*C_H) / (K_GAPDH_13BPG_i*K_GAPDH_NADH*K_GAPDH_Pi_i*K_GAPDH_13BPG_dash_i);

    r_GAPDH = (t0 - t1) / (t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11 + t12 + t13 + t14 + t15);
end
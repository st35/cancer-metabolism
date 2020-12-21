function r_LDH = Lactate_Dehydrogenase(C_PYR, C_NADH, C_LAC, C_NAD, Rate)
    V_LDH_mf = 8.66e3; V_LDH_mr = 2.17e3;
    K_LDH_PYR = 0.137; K_LDH_NADH = 7.43e-3; K_LDH_LAC = 1.07; K_LDH_NAD = 0.107;
    K_LDH_PYR_i = 0.228; K_LDH_NADH_i = 5.45e-3; K_LDH_LAC_i = 7.33; K_LDH_NAD_i = 0.503;
    K_LDH_PYR_i_dash = 0.101;

    V_LDH_mf = V_LDH_mf*1e2;
    V_LDH_mr = V_LDH_mr*1e2;

    V_LDH_mf = V_LDH_mf*Rate;
    V_LDH_mr = V_LDH_mr*Rate;

    t0 = ((V_LDH_mf*C_PYR*C_NADH) / (K_LDH_PYR*K_LDH_NADH_i)) - ((V_LDH_mr*C_LAC*C_NAD) / (K_LDH_LAC*K_LDH_NAD_i));
    t1 = (1.0 + ((K_LDH_NADH*C_PYR) / (K_LDH_NADH_i*K_LDH_PYR)) + ((K_LDH_NAD*C_LAC) / (K_LDH_NAD_i*K_LDH_LAC)))*(1.0 + (C_PYR / K_LDH_PYR_i_dash));
    t2 = (C_NADH / K_LDH_NADH_i) + (C_NAD / K_LDH_NAD_i);
    t3 = (C_NADH*C_PYR) / (K_LDH_NADH_i*K_LDH_PYR);
    t4 = (K_LDH_NAD*C_NADH*C_LAC) / (K_LDH_NAD_i*K_LDH_NADH_i*K_LDH_LAC);
    t5 = (K_LDH_NADH*C_NAD*C_PYR) / (K_LDH_NAD_i*K_LDH_NADH_i*K_LDH_PYR);
    t6 = (C_NAD*C_LAC) / (K_LDH_NAD_i*K_LDH_LAC);
    t7 = (C_NADH*C_PYR*C_LAC) / (K_LDH_NADH_i*K_LDH_PYR*K_LDH_LAC_i);
    t8 = (C_NAD*C_PYR*C_LAC) / (K_LDH_NAD_i*K_LDH_PYR_i*K_LDH_LAC);

    r_LDH = (t0) / (t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8);
end
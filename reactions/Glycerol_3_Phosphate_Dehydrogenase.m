function r_G3PD = Glycerol_3_Phosphate_Dehydrogenase(C_DHAP, C_NADH, C_G3P, C_NAD, C_ADP, C_ATP, C_F16BP, Rate)
    V_G3PD_mf = 67.0*60.0;
    K_G3PD_eq = 1e4;
    K_G3PD_NADH = 0.023; K_G3PD_DHAP = 0.54; K_G3PD_NAD = 0.93; K_G3PD_G3P = 1.2;
    K_G3PD_F16BP = 4.8; K_G3PD_ADP = 2.0; K_G3PD_ATP = 0.73;

    V_G3PD_mf = V_G3PD_mf*1.0e2;

    V_G3PD_mf = V_G3PD_mf*Rate;

    t0 = (V_G3PD_mf / (K_G3PD_NADH*K_G3PD_DHAP))*(C_NADH*C_DHAP - (C_NAD*C_G3P) / (K_G3PD_eq));
    t1 = (1.0 + (C_DHAP / K_G3PD_DHAP) + (C_G3P / K_G3PD_G3P));
    t2 = (1.0 + (C_NADH / K_G3PD_NADH) + (C_NAD / K_G3PD_NAD));
    i0 = (1.0 + (C_F16BP / K_G3PD_F16BP) + (C_ADP / K_G3PD_ADP) + (C_ATP / K_G3PD_ATP));

    r_G3PD = t0 / (t1*t2*i0);
end
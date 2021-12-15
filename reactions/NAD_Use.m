function r_NAD_Use = NAD_Use(C_NAD, C_NADH, Rate)
    V_NAD_Use = 0.426*3600.0;
    K_Use_NAD = 0.3*2.0; K_Use_NADH = 0.005;

    t0 = (V_NAD_Use*C_NAD / K_Use_NAD);
    t1 = (1.0 + (C_NAD / K_Use_NAD) + (C_NADH / K_Use_NADH));

    r_NAD_Use = t0 / t1;
end
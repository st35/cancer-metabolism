function r_ATP_Use = Cyto_ATP_Use(C_ATP, C_ADP, Rate)
    V_ATP_Use = 0.426*3600.0;
    k_ATP_Use = 3.0; k_ADP = 1.0e-1;

    V_ATP_Use = V_ATP_Use*1e3;

    V_ATP_Use = V_ATP_Use*Rate;

    t0 = (V_ATP_Use*C_ATP) / k_ATP_Use;
    t1 = 1.0 + (C_ATP / k_ATP_Use) + (C_ADP / k_ADP);

    r_ATP_Use = t0 / t1;
end
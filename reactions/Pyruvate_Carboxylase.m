function r_PC = Pyruvate_Carboxylase(M_PYR, M_ATP, M_CO2, M_OAA, M_ADP, Rate)
    V_PC_mf = 32.9*3600.0;
    K_PC_PYR = 0.49; K_PC_ATP = 0.13; K_PC_CO2 = 1.3;
    K_PC_OAA = 0.027; K_PC_ADP = 0.19;

    V_PC_mf = V_PC_mf*Rate;

    t0 = (V_PC_mf*M_PYR*M_ATP*M_CO2) / (K_PC_PYR*K_PC_ATP*K_PC_CO2);
    t1 = 1.0 + (M_PYR / K_PC_PYR) + (M_CO2 / K_PC_CO2) + (M_ATP / K_PC_ATP);
    t2 = ((M_PYR*M_CO2) / (K_PC_PYR*K_PC_CO2)) + ((M_PYR*M_ATP) / (K_PC_PYR*K_PC_ATP)) + ((M_CO2*M_ATP) / (K_PC_CO2*K_PC_ATP));
    t3 = ((M_PYR*M_CO2*M_ATP) / (K_PC_PYR*K_PC_CO2*K_PC_ATP));
    t4 = (M_OAA / K_PC_OAA) + (M_ADP / K_PC_ADP);

    r_PC = t0 / (t1 + t2 + t3 + t4);
end
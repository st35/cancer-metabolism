function r_SLC25A1 = Citrate_Malate_Antiport(M_CIT, M_MAL, C_CIT, C_MAL, Rate)
    V_SLC25A1_mf = 7.31e1*3600.0*1e3;

    V_SLC25A1_mf = V_SLC25A1_mf*Rate;

    r_SLC25A1 = V_SLC25A1_mf*(M_CIT*C_MAL - C_CIT*M_MAL);
end
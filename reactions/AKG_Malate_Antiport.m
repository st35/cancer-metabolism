function r_SLC25A11 = AKG_Malate_Antiport(M_AKG, M_MAL, C_AKG, C_MAL, Rate)
    V_SLC25A11_mf = 7.31e1*3600.0*1e3;

    V_SLC25A11_mf = V_SLC25A11_mf*Rate;

    r_SLC25A11 = V_SLC25A11_mf*(M_AKG*C_MAL - C_AKG*M_MAL);
end
function r_Phospho = Phospholipids_Synthesis(C_G3P, Rate)
    V_Phospho_mf = 67.0*60.0;
    K_Phospho_G3P = 1.2*2.0;

    V_Phospho_mf = V_Phospho_mf*Rate;

    r_Phospho = (V_Phospho_mf*C_G3P / K_Phospho_G3P) / (1.0 + (C_G3P / K_Phospho_G3P));
end
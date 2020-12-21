function r_ACC = Acetyl_CoA_Carboxylase(C_ACCOA, Rate)
    V_ACC_mf = 2.4e4;
    K_ACC_ACCOA = 28.8e-3;

    V_ACC_mf = V_ACC_mf*Rate;

    t0 = (V_ACC_mf*C_ACCOA / K_ACC_ACCOA);
    t1 = 1.0 + (C_ACCOA / K_ACC_ACCOA);

    r_ACC = t0 / t1;
end
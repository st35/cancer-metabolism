function r_IDH = Isocitrate_Dehydrogenase(M_NAD, M_ICIT, M_AKG, M_NADH, M_CO2, M_ADP, M_ATP, Rate)
    V_IDH_mf = 4.25e-1*1e3*3600.0; V_IDH_mr = V_IDH_mf / 1e2;
    n_H = 3.0; k_mB = 183.0e-3; k_mA = 74.0e-3; k_ib = 23.8e-3; k_iq = 29.0e-3;
    k_aADP = 50.0e-3; k_iATP = 91.0e-3;

    V_IDH_mf = V_IDH_mf*Rate;
    V_IDH_mr = V_IDH_mr*Rate;

    alpha_i = 1.0 + (k_aADP / M_ADP)*(1.0 + M_ATP / k_iATP);

    A = M_NAD; B = M_ICIT;
    P = M_AKG; Q = M_NADH; R = M_CO2;

    t0 = (V_IDH_mf*(A*B^n_H)) - (V_IDH_mr*(B^(n_H - 1))*P*Q*R);
    t1 = A*B^n_H + (k_mB^n_H)*alpha_i*A + k_mA*(B^n_H + k_ib^n_H*alpha_i + (alpha_i*(Q*B^n_H / k_iq)));
    
    r_IDH = t0 / t1;
end
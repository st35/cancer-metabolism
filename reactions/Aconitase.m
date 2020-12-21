function r_ACO = Aconitase(M_CIT, M_ICIT)
    V_ACO_mf = 3.21e-2*1e3*3600.0; V_ACO_mr = V_ACO_mf / 1e2;
    K_ACO_A = 1161.0e-3; K_ACO_B = 434.0e-3;

    t0 = ((V_ACO_mf*M_CIT) / K_ACO_A) - (V_ACO_mr*M_ICIT / K_ACO_B);
    t1 = 1.0 + (M_CIT / K_ACO_A) + (M_ICIT / K_ACO_B);

    r_ACO = t0 / t1;
end
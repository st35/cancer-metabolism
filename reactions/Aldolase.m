function r_ALD = Aldolase(C_F16BP, C_GAP, C_DHAP, C_23BPG)
    V_ALD_mf = 6.75e2; V_ALD_mr = 2.32e3;
    K_ALD_F16BP = 5.0e-2; K_ALD_GAP = 0.189; K_ALD_DHAP = 3.5e-2;
    K_ALD_23BPG = 1.5; K_ALD_F16BP_i = 1.98e-2; K_ALD_DHAP_i = 1.1e-2;

    V_ALD_mf = V_ALD_mf*1.0e2;
    V_ALD_mr = V_ALD_mr*1.0e2;
    
    t0 = ((V_ALD_mf*C_F16BP / K_ALD_F16BP)) - ((V_ALD_mr*C_GAP*C_DHAP) / (K_ALD_GAP*K_ALD_DHAP_i));
    t1 = 1.0 + (C_23BPG / K_ALD_23BPG) + (C_F16BP / K_ALD_F16BP) + (C_DHAP / K_ALD_DHAP_i);
    t2 = ((K_ALD_DHAP*C_F16BP*C_GAP) / (K_ALD_F16BP_i*K_ALD_GAP*K_ALD_DHAP_i));
    t3 = (C_GAP*C_DHAP) / (K_ALD_GAP*K_ALD_DHAP_i);
    t4 = ((K_ALD_DHAP*C_GAP) / (K_ALD_GAP*K_ALD_DHAP_i));
    t5 = 1.0 + (C_23BPG / K_ALD_23BPG);

    r_ALD = (t0) / (t1 + t2 + t3 + t4 + t5);
end
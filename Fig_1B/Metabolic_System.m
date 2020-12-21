function dxdt = Metabolic_System(type, N_Metabolites, N_Fluxes, Input, x)
    if length(x) ~= N_Metabolites
        'This is yet another Minute Maid fiasco. Damn it, Carl!'
    end

    C_MgATP = 2.69; C_MgADP = 0.46; C_G16BP = 0.0; C_23BPG = 3.0; C_GSH = 0.57; C_ATP = 0.31; C_Mg = 0.7; C_AMP = 0.03; C_Pi = 2.5; C_NAD = 0.299; C_NADH = 0.001;
    C_pH = 7.3; C_H = (10.0^(-C_pH))*1e3; C_ALA = 0.2; M_pH = 8.0; M_H = (10.0^(-M_pH))*1e3; C_NADP = 0.0014; C_NADPH = 0.0643; C_ADP = 0.54;
    M_COASH = 9.9e-2; M_NAD = (2.97e3)*(8.0 / 9.0); M_CO2 = 21.4; M_NADH = (2.97e3)*(1.0 / 9.0); M_OAA = 4.0e-2; M_ATP = 900.0e-4; M_ADP = 225.0e-4; M_AMP = 0.0;
    M_SCOA = 0.0; C_MAL = 1.0e-3; C_COA = 12.0e-2; M_Pi = 2.44; C_OAA = 1.0; M_GDP = 80.0e-3; M_GTP = 80.0e-5; M_COQ = 480.0e-2; M_COQH2 = 2.45e-4; M_NAD_Tot = 2.97e3;
    M_NH3 = 2.4e-3; C_CO2 = 21.4; C_ATP_Tot = 2.69 + 0.46 + 0.31 + 0.54;

    Rate_GLUT = Input(1);
    Rate_HK = 1.0;
    Rate_G6PD = 1.0;
    Rate_G3PD = 1.0;
    Rate_LDH = 1.0;
    Rate_PDH = 1.0;
    Rate_CIT = 1.0;
    Rate_IDH = 1.0;
    Rate_AKGD = 1.0;
    Rate_SCOAS = 1.0;
    Rate_SDH = 1.0;
    Rate_FUM = 1.0;
    Rate_MDH = 1.0;

    E_GLC = 1.0;

    C_GLC = x(1);
    C_G6P = x(2);
    C_F6P = x(3);
    C_6PGDL = 0.0;
    C_F16BP = x(4);
    C_DHAP = x(5);
    C_GAP = x(6);
    C_G3P = 0.0;
    C_13BPG = x(7);
    C_3PG = x(8);
    C_2PG = x(9);
    C_PEP = x(10);
    C_PYR = x(11);
    C_LAC = 0.0;
    M_PYR = x(12);
    M_ACCOA = x(13);
    M_CIT = x(14);
    M_ICIT = x(15);
    M_AKG = x(16);
    M_SCOA = x(17);
    M_SUC = x(18);
    M_FUM = x(19);
    M_MAL = x(20);
    M_OAA = x(21);
    C_CIT = 0.0;

    r_GLUT = Glucose_Transporter(E_GLC, C_GLC, Rate_GLUT);
    r_HK = Hexokinase(C_GLC, C_MgATP, C_G6P, C_MgADP, C_G16BP, C_23BPG, C_GSH, Rate_HK);
    r_GPI = Glucose_Phosphate_Isomerase(C_G6P, C_F6P);
    r_G6PD = Glucose_6_Phosphate_Dehydrogenase(C_G6P, C_NADP, C_6PGDL, C_NADPH, Rate_G6PD);
    r_PFK = Phosphofructokinase(C_F6P, C_MgATP, C_F16BP, C_MgADP, C_ATP, C_Mg, C_23BPG, C_AMP, C_G16BP, C_Pi, C_CIT);
    r_ALD = Aldolase(C_F16BP, C_GAP, C_DHAP, C_23BPG);
    r_TPI = Triose_Phosphate_Isomerase(C_DHAP, C_GAP);
    r_G3PD = Glycerol_3_Phosphate_Dehydrogenase(C_DHAP, C_NADH, C_G3P, C_NAD, C_ADP, C_ATP, C_F16BP, Rate_G3PD);
    r_GAPDH = Glyceraldehyde_3_Phosphate_Dehydrogenase(C_GAP, C_NAD, C_Pi, C_13BPG, C_NADH, C_H);
    r_PGK = Phosphoglycerate_Kinase(C_13BPG, C_MgADP, C_3PG, C_MgATP);
    r_PGM = Phosphoglycerate_Mutase(C_3PG, C_2PG);
    r_ENO = Enolase(C_2PG, C_Mg, C_PEP);
    r_PK = Pyruvate_Kinase(C_PEP, C_MgADP, C_PYR, C_MgATP, C_ATP, C_ALA, C_F16BP, C_G16BP);
    r_LDH = Lactate_Dehydrogenase(C_PYR, C_NADH, C_LAC, C_NAD, Rate_LDH);
    r_PIM = Mitochondrial_Pyruvate_Transport(C_PYR, C_H, M_PYR, M_H);
    r_PDH = Pyruvate_Dehydrogenase(M_PYR, M_COASH, M_NAD, M_CO2, M_ACCOA, M_NADH, Rate_PDH);
    r_CIT = Citrate_Synthase(M_ACCOA, M_OAA, M_CIT, M_COASH, M_ATP, M_ADP, M_AMP, M_SCOA, Rate_CIT);
    r_ACO = Aconitase(M_CIT, M_ICIT);
    r_IDH = Isocitrate_Dehydrogenase(M_NAD, M_ICIT, M_AKG, M_NADH, M_CO2, M_ADP, M_ATP, Rate_IDH);
    r_AKGD = AKG_Dehydrogenase(M_AKG, M_COASH, M_NAD, M_SCOA, M_NADH, Rate_AKGD);
    r_SCOAS = Succinyl_CoA_Synthetase(M_GDP, M_SCOA, M_Pi, M_COASH, M_SUC, M_GTP, Rate_SCOAS);
    r_SDH = Succinate_Dehydrogenase(M_SUC, M_COQ, M_FUM, M_COQH2, M_OAA, Rate_SDH);
    r_FUM = Fumarase(M_FUM, M_MAL, M_CIT, M_ATP, M_ADP, M_GTP, M_GDP, Rate_FUM);
    r_MDH = Malate_Dehydrogenase(M_MAL, M_NAD, M_OAA, M_NADH, Rate_MDH);

    if type == 0
        dxdt = zeros(N_Metabolites, 1);
        dxdt(1) = r_GLUT - r_HK;
        dxdt(2) = r_HK - r_GPI - r_G6PD;
        dxdt(3) = r_GPI - r_PFK;
        dxdt(4) = r_PFK - r_ALD;
        dxdt(5) = r_ALD - r_TPI - r_G3PD;
        dxdt(6) = r_ALD + r_TPI - r_GAPDH;
        dxdt(7) = r_GAPDH - r_PGK;
        dxdt(8) = r_PGK - r_PGM;
        dxdt(9) = r_PGM - r_ENO;
        dxdt(10) = r_ENO - r_PK;
        dxdt(11) = r_PK - r_LDH - r_PIM;
        dxdt(12) = r_PIM - r_PDH;
        dxdt(13) = r_PDH - r_CIT;
        dxdt(14) = r_CIT - r_ACO;
        dxdt(15) = r_ACO - r_IDH;
        dxdt(16) = r_IDH - r_AKGD;
        dxdt(17) = r_AKGD - r_SCOAS;
        dxdt(18) = r_SCOAS - r_SDH;
        dxdt(19) = r_SDH - r_FUM;
        dxdt(20) = r_FUM - r_MDH;
        dxdt(21) = r_MDH - r_CIT;
    else
        dxdt = zeros(N_Fluxes, 1);
        dxdt(1) = r_GLUT;
        dxdt(2) = r_HK;
        dxdt(3) = r_GPI;
        dxdt(4) = r_G6PD;
        dxdt(5) = r_PFK;
        dxdt(6) = r_ALD;
        dxdt(7) = r_TPI;
        dxdt(8) = r_G3PD;
        dxdt(9) = r_GAPDH;
        dxdt(10) = r_PGK;
        dxdt(11) = r_PGM;
        dxdt(12) = r_ENO;
        dxdt(13) = r_PK;
        dxdt(14) = r_LDH;
        dxdt(15) = r_PIM;
        dxdt(16) = r_PDH;
        dxdt(17) = r_CIT;
        dxdt(18) = r_ACO;
        dxdt(19) = r_IDH;
        dxdt(20) = r_AKGD;
        dxdt(21) = r_SCOAS;
        dxdt(22) = r_SDH;
        dxdt(23) = r_FUM;
        dxdt(24) = r_MDH;
    end
end
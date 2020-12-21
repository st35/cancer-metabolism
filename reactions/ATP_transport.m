function r_ANT = ATP_transport(M_ATP, Rate)
    k_ANT = 0.426*1e3*3600.0;

    k_ANT = k_ANT*Rate;

    r_ANT = k_ANT*M_ATP;
end
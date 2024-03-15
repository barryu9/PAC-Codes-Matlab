function C = pac_fano_updatePSBack(pac_params, i_minus1, s_max, u_esti, C)
    k = 2^s_max;
    for i = i_minus1 + 1 - k:i_minus1
        C = update_C(pac_params, i, C, u_esti(i+1));
    end
end
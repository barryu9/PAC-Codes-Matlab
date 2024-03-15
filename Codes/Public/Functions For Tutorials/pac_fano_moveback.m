function [Threshold, jj, visited_before, P, C] = pac_fano_moveback(pac_params,rp, bmetric, bmetric_cut, j, Threshold, delta, u_esti, P, C,llr, Delta)
info_bits_indices=rp.info_bits_indices;
while 1
    follow_other_branch = false;
    if j >= 1
        mu_pre = bmetric(j);
    else
        mu_pre = 0;
    end

    jj = j;
    for k = j-1 :-1:0
        mu_pre = bmetric(k+1);
        if mu_pre >= Threshold
            if bmetric_cut(k+1) >= Threshold
                jj = k;
                tmp = bmetric(k+1);
                bmetric(k+1) = bmetric_cut(k+1);
                bmetric_cut(k+1) = tmp;
                if delta(k+1) == 0
                    follow_other_branch = true;
                    break;
                elseif k == 0
                    follow_other_branch = true;
                    jj = 0;
                    break;
                end
            end
        end
        if k == 0
            mu_pre = -inf;
        end
    end

    i_cur = info_bits_indices(j+1) - 1;
    if (mu_pre >= Threshold) && (jj ~= 0 || follow_other_branch == true)
        i_start = info_bits_indices(jj+1) - 1;
        [P, C] = pac_fano_updateLLRsPSs(pac_params,llr, i_start, i_cur, u_esti, P, C);

        if delta(jj+1) == 0
            visited_before = true;
            return;
        elseif jj == 0
            Threshold = Threshold - Delta;
            visited_before = false;
            return;
        end
    else
        Threshold = Threshold - Delta;
        visited_before = false;
        jj = j;
        return;
    end
end

end
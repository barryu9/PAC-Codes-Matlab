function [P, C] = pac_fano_updateLLRsPSs(pac_params,llr, i_start, i_cur, u_esti, P, C)
n = pac_params.n;
if mod(i_cur, 2) ~= 0
    i_cur = i_cur - 1;
end
if mod(i_start, 2) ~= 0
    i_start = i_start - 1;
end
s_start = ffs(i_start, n);
s_max = smax(i_start, i_cur, n);

if s_start <= s_max
    i_minus1 = find_sMaxPos(s_start, s_max, i_start, n);
    for i = i_minus1:i_start
        P = update_P(pac_params, llr,i, P, C);
        C = update_C(pac_params, i, C, u_esti(i+1));
    end
else
    P = update_P(pac_params, i_start, P, C);
end
end
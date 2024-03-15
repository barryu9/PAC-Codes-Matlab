function [d_esti] = pac_fano_decoder_slow(pac_params, rp,llr, Delta)
% Init;
N = pac_params.N;
k = pac_params.k;
conv_depth=pac_params.conv_depth;
gen = pac_params.gen;
pe = zeros(1,N);
% pe=rp.pe;
c_state = zeros(conv_depth-1, 1);
curr_state = zeros(conv_depth-1, k); %0~k-1,1:|g|-1;
phi = 0;
psi = 0;
Threshold = 0;

P = zeros(N-1, 1);
C = zeros(N-1, 2);
B = sum(log2(1-pe));
alphaq = 1;
info_bits_indices = rp.info_bits_indices;
frozen_bits_mask = rp.frozen_bits_mask;
delta = zeros(k, 1);
mu = zeros(N, 1);
bmetric = zeros(k, 1);
bmetric_cut = zeros(k, 1);
u_esti = zeros(N, 1);
v_esti = zeros(N, 1);
visited_before = false;
%while not end of tree do
while phi < N
    P = update_P(pac_params, llr,phi, P, C);
    if frozen_bits_mask(phi+1) == 1
        [u_esti(phi+1), c_state] = conv1bTrans(0, c_state, gen);
        if phi == 0
            mu(phi+1) = B + m_func(P(1), u_esti(phi+1)) - alphaq * log2(1-pe(phi+1));
        else
            mu(phi+1) = mu(phi) + m_func(P(1), u_esti(phi+1)) - alphaq * log2(1-pe(phi+1));
        end
        curr_state(:, psi+1) = c_state;
        C = update_C(pac_params, phi, C, u_esti(phi+1));
        phi = phi + 1;
    else

        % look forward to best node
        [u_left, c_state_left] = conv1bTrans(0, c_state, gen);
        [u_right, c_state_right] = conv1bTrans(1, c_state, gen);
        if phi == 0
            mu_left = B + m_func(P(1), u_left) - alphaq * log2(1-pe(phi+1));
            mu_right = B + m_func(P(1), u_right) - alphaq * log2(1-pe(phi+1));
        else
            mu_left = mu(phi) + m_func(P(1), u_left) - alphaq * log2(1-pe(phi+1));
            mu_right = mu(phi) + m_func(P(1), u_right) - alphaq * log2(1-pe(phi+1));
        end
        if mu_left > mu_right
            mu_max = mu_left;
            mu_min = mu_right;

        else
            mu_max = mu_right;
            mu_min = mu_left;
        end

        if mu_max >= Threshold
            if visited_before == false
                if mu_left > mu_right
                    v_esti(phi+1) = 0;
                    u_esti(phi+1) = u_left;
                else
                    v_esti(phi+1) = 1;
                    u_esti(phi+1) = u_right;
                end
                bmetric(psi+1) = mu_max;
                bmetric_cut(psi+1) = mu_min;
                mu(phi+1) = mu_max;
                delta(psi+1) = 0;
                C = update_C(pac_params, phi, C, u_esti(phi+1));
                curr_state(:, psi+1) = c_state;
                if (v_esti(phi+1) == 0)
                    c_state = c_state_left;
                else
                    c_state = c_state_right;
                end
                phi = phi + 1;
                psi = psi + 1;
            else
                if mu_min > Threshold
                    if mu_left < mu_right
                        v_esti(phi+1) = 0;
                        u_esti(phi+1) = u_left;
                    else
                        v_esti(phi+1) = 1;
                        u_esti(phi+1) = u_right;
                    end
                    bmetric(psi+1) = mu_min;
                    bmetric_cut(psi+1) = mu_max;
                    mu(phi+1) = mu_min;
                    C = update_C(pac_params, phi, C, u_esti(phi+1));
                    curr_state(:, psi+1) = c_state;
                    if (v_esti(phi+1) == 0)
                        c_state = c_state_left;
                    else
                        c_state = c_state_right;
                    end
                    delta(psi+1) = 1;
                    phi = phi + 1;
                    psi = psi + 1;
                    visited_before = false;
                else
                    if psi == 0
                        Threshold = Threshold - Delta;
                        visited_before = false;
                    else
                        curr_state(:, psi+1) = c_state;
                        [Threshold, psi, visited_before, P, C] = pac_fano_moveback(pac_params, bmetric, bmetric_cut, psi, Threshold, delta, u_esti, P, C, llr, Delta);
                        if psi == 0
                            c_state = zeros(conv_depth-1, 1);
                        else
                            c_state = curr_state(:, psi+1);
                        end
                        phi = info_bits_indices(psi+1) - 1;
                    end
                end
            end
        else
            if psi == 0
                while mu_max < Threshold
                    Threshold = Threshold - Delta;
                end
                visited_before = false;
            else
                curr_state(:, psi+1) = c_state;
                [Threshold, psi, visited_before, P, C] = pac_fano_moveback(pac_params,rp, bmetric, bmetric_cut, psi, Threshold, delta, u_esti, P, C, llr, Delta);
                if psi == 0
                    c_state = zeros(conv_depth-1, 1);
                else
                    c_state = curr_state(:, psi+1);
                end
                phi = info_bits_indices(psi+1) - 1;
            end
        end
    end
end
d_esti = v_esti(info_bits_indices);
end
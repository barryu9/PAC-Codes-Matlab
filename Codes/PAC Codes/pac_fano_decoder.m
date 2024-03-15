function [d_esti] = pac_fano_decoder(pac_params, rp,llr_batch, Delta)
% Init;
N = pac_params.N;
k = pac_params.k;
conv_depth=pac_params.conv_depth;
gen = pac_params.gen;

llr_layer_vec = pac_params.llr_layer_vec;
lambda_offset= pac_params.lambda_offset;
bit_layer_vec= pac_params.bit_layer_vec;
n=pac_params.n;

info_bits_indices = rp.info_bits_indices;
frozen_bits_mask = rp.frozen_bits_mask;
% pe = rp.pe;
pe = zeros(1,N);

batchsize = size(llr_batch,2);
d_esti = zeros(k,batchsize);
alphaq = 1;

for batch_index = 1:batchsize

    c_state = zeros(conv_depth-1, 1);
    curr_state = zeros(conv_depth-1, k); %0~k-1,1:|g|-1;
    phi = 0;
    psi = 0;
    Threshold = 0;
    P = zeros(N-1, 1);
    C = zeros(N-1, 2);
    B = sum(alphaq*log2(1-pe));
    delta = zeros(k, 1);
    mu = zeros(N, 1);
    bmetric = zeros(k, 1);
    bmetric_cut = zeros(k, 1);
    u_esti = zeros(N, 1);
    v_esti = zeros(N, 1);
    visited_before = false;
    llr = llr_batch(:,batch_index);

    %while not end of tree do
    while phi < N
        %     P = update_P(pac_params, llr,phi, P, C);
        % -------------------------------------------------
        layer = llr_layer_vec(phi + 1);
        switch phi %Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits.
            case 0
                index_1 = lambda_offset(n);
                for beta = 0:index_1 - 1
                    P(beta + index_1) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                    %             operation_count_P=operation_count_P+1;
                end
                for i_layer = n - 2:-1:0
                    index_1 = lambda_offset(i_layer+1);
                    index_2 = lambda_offset(i_layer+2);
                    for beta = 0:index_1 - 1
                        P(beta + index_1) = sign(P(beta + index_2)) * ...
                            sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                        %                 operation_count_P=operation_count_P+1;
                    end
                end
            case N / 2
                index_1 = lambda_offset(n);
                for beta = 0:index_1 - 1
                    x_tmp = C(beta + index_1,1);
                    P(beta + index_1) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                    %             operation_count_P=operation_count_P+1;

                end
                for i_layer = n - 2:-1:0
                    index_1 = lambda_offset(i_layer+1);
                    index_2 = lambda_offset(i_layer+2);
                    for beta = 0:index_1 - 1
                        P(beta + index_1) = sign(P(beta + index_2)) * ...
                            sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                        %                 operation_count_P=operation_count_P+1;

                    end
                end
            otherwise
                index_1 = lambda_offset(layer+1);
                index_2 = lambda_offset(layer+2);
                for beta = 0:index_1 - 1
                    P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * P(beta + index_2) + ...
                        P(beta + index_1 + index_2);
                    %             operation_count_P=operation_count_P+1;

                end
                for i_layer = layer - 1:-1:0
                    index_1 = lambda_offset(i_layer+1);
                    index_2 = lambda_offset(i_layer+2);
                    for beta = 0:index_1 - 1
                        P(beta + index_1) = sign(P(beta + index_2)) * ...
                            sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), ...
                            abs(P(beta + index_1 + index_2)));
                        %                 operation_count_P=operation_count_P+1;

                    end
                end
        end
        %   -----------------------------------------------------------------------
        if frozen_bits_mask(phi+1) == 1
            %   [u_esti(phi+1), c_state] = conv1bTrans(0, c_state, gen);
            %   -----------------------------------------------------------------------
            u_esti(phi+1)=mod(0*gen(1),2);
            for j=2:length(gen)
                if(gen(j)==1)
                    u_esti(phi+1)=mod(u_esti(phi+1)+c_state(j-1),2);
                end
            end
            for i=conv_depth-1:-1:2
                c_state(i)=c_state(i-1);
            end
            c_state(1)=0;
            %   -----------------------------------------------------------------------
            if phi == 0
                mu(phi+1) = 0  - log2(1+exp(-(1 - 2 * u_esti(phi+1))*P(1))) - alphaq * log2(1-pe(phi+1));
            else
                mu(phi+1) = mu(phi) - log2(1+exp(-(1 - 2 * u_esti(phi+1))*P(1))) - alphaq * log2(1-pe(phi+1));
            end
            curr_state(:, psi+1) = c_state;
            %         C = update_C(pac_params, phi, C, u_esti(phi+1));
            %   -----------------------------------------------------------------------
            phi_mod_2 = mod(phi, 2);
            C(1, 1+phi_mod_2) = u_esti(phi+1);

            if (phi_mod_2 == 1) && (phi ~= N - 1)
                layer = bit_layer_vec(phi + 1);
                for i_layer = 0:layer - 1
                    index_1 = lambda_offset(i_layer+1);
                    index_2 = lambda_offset(i_layer+2);
                    for beta = index_1:2 * index_1 - 1
                        C(beta + index_1, 2) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                        C(beta + index_2, 2) = C(beta, 2);
                        %             operation_count_C=operation_count_C+1;

                    end
                end
                index_1 = lambda_offset(layer+1);
                index_2 = lambda_offset(layer+2);
                for beta = index_1:2 * index_1 - 1
                    C(beta + index_1, 1) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                    C(beta + index_2, 1) = C(beta, 2);
                    %         operation_count_C=operation_count_C+1;
                end
            end
            %   -----------------------------------------------------------------------
            phi = phi + 1;
        else

            % look forward to best node
            %         [u_left, c_state_left] = conv1bTrans(0, c_state, gen);
            u_left=mod(0*gen(1),2);
            for j=2:length(gen)
                if(gen(j)==1)
                    u_left=mod(u_left+c_state(j-1),2);
                end
            end
            for i=conv_depth-1:-1:2
                c_state_left(i)=c_state(i-1);
            end
            c_state_left(1)=0;

            %         [u_right, c_state_right] = conv1bTrans(1, c_state, gen);
            u_right=mod(1*gen(1),2);
            for j=2:length(gen)
                if(gen(j)==1)
                    u_right=mod(u_right+c_state(j-1),2);
                end
            end
            for i=conv_depth-1:-1:2
                c_state_right(i)=c_state(i-1);
            end
            c_state_right(1)=1;

            if phi == 0
                mu_left = B -log2(1+exp(-(1 - 2 * u_left)*P(1))) - alphaq * log2(1-pe(phi+1));
                mu_right = B -log2(1+exp(-(1 - 2 * u_right)*P(1)))  - alphaq * log2(1-pe(phi+1));
            else
                mu_left = mu(phi) -log2(1+exp(-(1 - 2 * u_left)*P(1))) - alphaq * log2(1-pe(phi+1));
                mu_right = mu(phi) -log2(1+exp(-(1 - 2 * u_right)*P(1))) - alphaq * log2(1-pe(phi+1));
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
                    %   C = update_C(pac_params, phi, C, u_esti(phi+1));
                    %   -----------------------------------------------------------------------
                    phi_mod_2 = mod(phi, 2);
                    C(1, 1+phi_mod_2) = u_esti(phi+1);

                    if (phi_mod_2 == 1) && (phi ~= N - 1)
                        layer = bit_layer_vec(phi + 1);
                        for i_layer = 0:layer - 1
                            index_1 = lambda_offset(i_layer+1);
                            index_2 = lambda_offset(i_layer+2);
                            for beta = index_1:2 * index_1 - 1
                                C(beta + index_1, 2) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                                C(beta + index_2, 2) = C(beta, 2);
                                %             operation_count_C=operation_count_C+1;

                            end
                        end
                        index_1 = lambda_offset(layer+1);
                        index_2 = lambda_offset(layer+2);
                        for beta = index_1:2 * index_1 - 1
                            C(beta + index_1, 1) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                            C(beta + index_2, 1) = C(beta, 2);
                            %         operation_count_C=operation_count_C+1;
                        end
                    end
                    %   -----------------------------------------------------------------------
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
                        %   C = update_C(pac_params, phi, C, u_esti(phi+1));
                        %   -----------------------------------------------------------------------
                        phi_mod_2 = mod(phi, 2);
                        C(1, 1+phi_mod_2) = u_esti(phi+1);

                        if (phi_mod_2 == 1) && (phi ~= N - 1)
                            layer = bit_layer_vec(phi + 1);
                            for i_layer = 0:layer - 1
                                index_1 = lambda_offset(i_layer+1);
                                index_2 = lambda_offset(i_layer+2);
                                for beta = index_1:2 * index_1 - 1
                                    C(beta + index_1, 2) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                                    C(beta + index_2, 2) = C(beta, 2);
                                    %             operation_count_C=operation_count_C+1;

                                end
                            end
                            index_1 = lambda_offset(layer+1);
                            index_2 = lambda_offset(layer+2);
                            for beta = index_1:2 * index_1 - 1
                                C(beta + index_1, 1) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                                C(beta + index_2, 1) = C(beta, 2);
                                %         operation_count_C=operation_count_C+1;
                            end
                        end
                        %   -----------------------------------------------------------------------
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
                            %                         [Threshold, psi, visited_before, P, C] = pac_fano_moveback(pac_params,rp, bmetric, bmetric_cut, psi, Threshold, delta, u_esti, P, C, llr, Delta);
                            %   -----------------------------------------------------------------------
                            while 1
                                follow_other_branch = false;
                                if psi >= 1
                                    mu_pre = bmetric(psi);
                                else
                                    mu_pre = 0;
                                end

                                psi_temp = psi;
                                for kk = psi - 1:-1:0
                                    mu_pre = bmetric(kk+1);
                                    if mu_pre >= Threshold
                                        if bmetric_cut(kk+1) >= Threshold
                                            psi_temp = kk;
                                            tmp = bmetric(kk+1);
                                            bmetric(kk+1) = bmetric_cut(kk+1);
                                            bmetric_cut(k+1) = tmp;
                                            if delta(kk+1) == 0
                                                follow_other_branch = true;
                                                break;
                                            elseif kk == 0
                                                follow_other_branch = true;
                                                psi_temp = 0;
                                                break;
                                            end
                                        end
                                    end
                                    if kk == 0
                                        mu_pre = -100;
                                    end
                                end

                                i_cur = info_bits_indices(psi+1) - 1;
                                if (mu_pre >= Threshold) && (psi_temp ~= 0 || follow_other_branch == true)
                                    i_start = info_bits_indices(psi_temp+1) - 1;
                                    %[P, C] = pac_fano_updateLLRsPSs(pac_params,llr, i_start, i_cur, u_esti, P, C);
                                    %-----------------------------------------------------------
                                    if mod(i_cur, 2) ~= 0
                                        i_cur = i_cur - 1;
                                    end
                                    if mod(i_start, 2) ~= 0
                                        i_start = i_start - 1;
                                    end
                                    %s_start = ffs(i_start, n);
                                    %-----------------------------------------------------------
                                    if i_start == 0
                                        s_start = n - 1;
                                        break;
                                    end
                                    bin = abs(dec2bin(i_start, n)) - 48;
                                    s_start = 0;
                                    for j = length(bin):-1:1
                                        if (bin(j) == 1)
                                            break;
                                        end
                                        s_start = s_start + 1;
                                    end
                                    %-----------------------------------------------------------
                                    %s_max = smax(i_start, i_cur, n);
                                    %-----------------------------------------------------------
                                    s_max = 0;
                                    for im = i_start:i_cur
                                        %temp = ffs(im, n);
                                        %-----------------------------------------------------------
                                        if im == 0
                                            temp = n - 1;
                                            break;
                                        end
                                        bin = abs(dec2bin(im, n)) - 48;
                                        temp = 0;
                                        for j = length(bin):-1:1
                                            if (bin(j) == 1)
                                                break;
                                            end
                                            temp = temp + 1;
                                        end
                                        %-----------------------------------------------------------
                                        if (temp > s_max)
                                            s_max = temp;
                                        end
                                    end
                                    %-----------------------------------------------------------

                                    if s_start <= s_max
                                        %i_minus1 = find_sMaxPos(s_start, s_max, i_start, n);
                                        %--------------------------------------------------------
                                        i_minus1 = i_start;
                                        ss = s_start;
                                        while ss < s_max
                                            i_minus1 = i_minus1 - 2;
                                            if i_minus1 > 0
                                                %ss = ffs(i_minus1, n);
                                                %-----------------------------------------------
                                                if i_minus1 == 0
                                                    ss = n - 1;
                                                    break;
                                                end
                                                bin = abs(dec2bin(i_minus1, n)) - 48;
                                                ss = 0;
                                                for j = length(bin):-1:1
                                                    if (bin(j) == 1)
                                                        break;
                                                    end
                                                    ss = ss + 1;
                                                end
                                                %-----------------------------------------------
                                            else
                                                ss = n;
                                            end
                                        end
                                        %-------------------------------------------------------
                                        for i = i_minus1:i_start
                                            %P = update_P(pac_params, llr,i, P, C);
                                            %-------------------------------------------------------
                                            layer = llr_layer_vec(i + 1);
                                            switch i %Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits.
                                                case 0
                                                    index_1 = lambda_offset(n);
                                                    for beta = 0:index_1 - 1
                                                        P(beta + index_1) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                                                        %             operation_count_P=operation_count_P+1;
                                                    end
                                                    for i_layer = n - 2:-1:0
                                                        index_1 = lambda_offset(i_layer+1);
                                                        index_2 = lambda_offset(i_layer+2);
                                                        for beta = 0:index_1 - 1
                                                            P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                                sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                            %                 operation_count_P=operation_count_P+1;
                                                        end
                                                    end
                                                case N / 2
                                                    index_1 = lambda_offset(n);
                                                    for beta = 0:index_1 - 1
                                                        x_tmp = C(beta + index_1,1);
                                                        P(beta + index_1) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                                                        %             operation_count_P=operation_count_P+1;

                                                    end
                                                    for i_layer = n - 2:-1:0
                                                        index_1 = lambda_offset(i_layer+1);
                                                        index_2 = lambda_offset(i_layer+2);
                                                        for beta = 0:index_1 - 1
                                                            P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                                sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                            %                 operation_count_P=operation_count_P+1;

                                                        end
                                                    end
                                                otherwise
                                                    index_1 = lambda_offset(layer+1);
                                                    index_2 = lambda_offset(layer+2);
                                                    for beta = 0:index_1 - 1
                                                        P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * P(beta + index_2) + ...
                                                            P(beta + index_1 + index_2);
                                                        %             operation_count_P=operation_count_P+1;

                                                    end
                                                    for i_layer = layer - 1:-1:0
                                                        index_1 = lambda_offset(i_layer+1);
                                                        index_2 = lambda_offset(i_layer+2);
                                                        for beta = 0:index_1 - 1
                                                            P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                                sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), ...
                                                                abs(P(beta + index_1 + index_2)));
                                                            %                 operation_count_P=operation_count_P+1;

                                                        end
                                                    end
                                            end
                                            %-------------------------------------------------------
                                            %C = update_C(pac_params, i, C, u_esti(i+1));
                                            %-------------------------------------------------------
                                            phi_mod_2 = mod(i, 2);
                                            C(1, 1+phi_mod_2) = u_esti(i+1);

                                            if (phi_mod_2 == 1) && (i ~= N - 1)
                                                layer = bit_layer_vec(i + 1);
                                                for i_layer = 0:layer - 1
                                                    index_1 = lambda_offset(i_layer+1);
                                                    index_2 = lambda_offset(i_layer+2);
                                                    for beta = index_1:2 * index_1 - 1
                                                        C(beta + index_1, 2) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                                                        C(beta + index_2, 2) = C(beta, 2);
                                                        %             operation_count_C=operation_count_C+1;

                                                    end
                                                end
                                                index_1 = lambda_offset(layer+1);
                                                index_2 = lambda_offset(layer+2);
                                                for beta = index_1:2 * index_1 - 1
                                                    C(beta + index_1, 1) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                                                    C(beta + index_2, 1) = C(beta, 2);
                                                    %         operation_count_C=operation_count_C+1;
                                                end
                                            end
                                            %-------------------------------------------------------

                                        end
                                    else
                                        %P = update_P(pac_params, i_start, P, C);
                                        %-------------------------------------------------------
                                        layer = llr_layer_vec(i_start + 1);
                                        switch i_start %Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits.
                                            case 0
                                                index_1 = lambda_offset(n);
                                                for beta = 0:index_1 - 1
                                                    P(beta + index_1) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                                                    %             operation_count_P=operation_count_P+1;
                                                end
                                                for i_layer = n - 2:-1:0
                                                    index_1 = lambda_offset(i_layer+1);
                                                    index_2 = lambda_offset(i_layer+2);
                                                    for beta = 0:index_1 - 1
                                                        P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                            sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                        %                 operation_count_P=operation_count_P+1;
                                                    end
                                                end
                                            case N / 2
                                                index_1 = lambda_offset(n);
                                                for beta = 0:index_1 - 1
                                                    x_tmp = C(beta + index_1,1);
                                                    P(beta + index_1) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                                                    %             operation_count_P=operation_count_P+1;

                                                end
                                                for i_layer = n - 2:-1:0
                                                    index_1 = lambda_offset(i_layer+1);
                                                    index_2 = lambda_offset(i_layer+2);
                                                    for beta = 0:index_1 - 1
                                                        P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                            sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                        %                 operation_count_P=operation_count_P+1;

                                                    end
                                                end
                                            otherwise
                                                index_1 = lambda_offset(layer+1);
                                                index_2 = lambda_offset(layer+2);
                                                for beta = 0:index_1 - 1
                                                    P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * P(beta + index_2) + ...
                                                        P(beta + index_1 + index_2);
                                                    %             operation_count_P=operation_count_P+1;

                                                end
                                                for i_layer = layer - 1:-1:0
                                                    index_1 = lambda_offset(i_layer+1);
                                                    index_2 = lambda_offset(i_layer+2);
                                                    for beta = 0:index_1 - 1
                                                        P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                            sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), ...
                                                            abs(P(beta + index_1 + index_2)));
                                                        %                 operation_count_P=operation_count_P+1;

                                                    end
                                                end
                                        end
                                        %-------------------------------------------------------
                                    end
                                    %-----------------------------------------------------------

                                    if delta(psi_temp+1) == 0
                                        visited_before = true;
                                        break;
                                    elseif psi_temp == 0
                                        Threshold = Threshold - Delta;
                                        visited_before = false;
                                        break;

                                    end
                                else
                                    Threshold = Threshold - Delta;
                                    visited_before = false;
                                    break;
                                end
                            end
                            psi = psi_temp;

                            %   -----------------------------------------------------------------------

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
                    %                 [Threshold, psi, visited_before, P, C] = pac_fano_moveback(pac_params,rp, bmetric, bmetric_cut, psi, Threshold, delta, u_esti, P, C, llr, Delta);
                    %------------------------------------------------------------------------------
                    while 1
                        follow_other_branch = false;
                        if psi >= 1
                            mu_pre = bmetric(psi);
                        else
                            mu_pre = 0;
                        end

                        psi_temp = psi;
                        for kk = psi - 1:-1:0
                            mu_pre = bmetric(kk+1);
                            if mu_pre >= Threshold
                                if bmetric_cut(kk+1) >= Threshold
                                    psi_temp = kk;
                                    tmp = bmetric(kk+1);
                                    bmetric(kk+1) = bmetric_cut(kk+1);
                                    bmetric_cut(kk+1) = tmp;
                                    if delta(kk+1) == 0
                                        follow_other_branch = true;
                                        break;
                                    elseif kk == 0
                                        follow_other_branch = true;
                                        psi_temp = 0;
                                        break;
                                    end
                                end
                            end
                            if kk == 0
                                mu_pre = -100;
                            end
                        end

                        i_cur = info_bits_indices(psi+1) - 1;
                        if (mu_pre >= Threshold) && (psi_temp ~= 0 || follow_other_branch == true)
                            i_start = info_bits_indices(psi_temp+1) - 1;
                            %[P, C] = pac_fano_updateLLRsPSs(pac_params,llr, i_start, i_cur, u_esti, P, C);
                            %-----------------------------------------------------------
                            if mod(i_cur, 2) ~= 0
                                i_cur = i_cur - 1;
                            end
                            if mod(i_start, 2) ~= 0
                                i_start = i_start - 1;
                            end
                            %s_start = ffs(i_start, n);
                            %-----------------------------------------------------------
                            if i_start == 0
                                s_start = n - 1;
                                break;
                            end
                            bin = abs(dec2bin(i_start, n)) - 48;
                            s_start = 0;
                            for j = length(bin):-1:1
                                if (bin(j) == 1)
                                    break;
                                end
                                s_start = s_start + 1;
                            end
                            %-----------------------------------------------------------
                            %s_max = smax(i_start, i_cur, n);
                            %-----------------------------------------------------------
                            s_max = 0;
                            for im = i_start:i_cur
                                %temp = ffs(im, n);
                                %-----------------------------------------------------------
                                if im == 0
                                    temp = n - 1;
                                    break;
                                end
                                bin = abs(dec2bin(im, n)) - 48;
                                temp = 0;
                                for j = length(bin):-1:1
                                    if (bin(j) == 1)
                                        break;
                                    end
                                    temp = temp + 1;
                                end
                                %-----------------------------------------------------------
                                if (temp > s_max)
                                    s_max = temp;
                                end
                            end
                            %-----------------------------------------------------------

                            if s_start <= s_max
                                %i_minus1 = find_sMaxPos(s_start, s_max, i_start, n);
                                %--------------------------------------------------------
                                i_minus1 = i_start;
                                ss = s_start;
                                while ss < s_max
                                    i_minus1 = i_minus1 - 2;
                                    if i_minus1 > 0
                                        %ss = ffs(i_minus1, n);
                                        %-----------------------------------------------
                                        if i_minus1 == 0
                                            ss = n - 1;
                                            break;
                                        end
                                        bin = abs(dec2bin(i_minus1, n)) - 48;
                                        ss = 0;
                                        for j = length(bin):-1:1
                                            if (bin(j) == 1)
                                                break;
                                            end
                                            ss = ss + 1;
                                        end
                                        %-----------------------------------------------
                                    else
                                        ss = n;
                                    end
                                end
                                %-------------------------------------------------------
                                for i = i_minus1:i_start
                                    %P = update_P(pac_params, llr,i, P, C);
                                    %-------------------------------------------------------
                                    layer = llr_layer_vec(i + 1);
                                    switch i %Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits.
                                        case 0
                                            index_1 = lambda_offset(n);
                                            for beta = 0:index_1 - 1
                                                P(beta + index_1) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                                                %             operation_count_P=operation_count_P+1;
                                            end
                                            for i_layer = n - 2:-1:0
                                                index_1 = lambda_offset(i_layer+1);
                                                index_2 = lambda_offset(i_layer+2);
                                                for beta = 0:index_1 - 1
                                                    P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                        sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                    %                 operation_count_P=operation_count_P+1;
                                                end
                                            end
                                        case N / 2
                                            index_1 = lambda_offset(n);
                                            for beta = 0:index_1 - 1
                                                x_tmp = C(beta + index_1,1);
                                                P(beta + index_1) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                                                %             operation_count_P=operation_count_P+1;

                                            end
                                            for i_layer = n - 2:-1:0
                                                index_1 = lambda_offset(i_layer+1);
                                                index_2 = lambda_offset(i_layer+2);
                                                for beta = 0:index_1 - 1
                                                    P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                        sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                    %                 operation_count_P=operation_count_P+1;

                                                end
                                            end
                                        otherwise
                                            index_1 = lambda_offset(layer+1);
                                            index_2 = lambda_offset(layer+2);
                                            for beta = 0:index_1 - 1
                                                P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * P(beta + index_2) + ...
                                                    P(beta + index_1 + index_2);
                                                %             operation_count_P=operation_count_P+1;

                                            end
                                            for i_layer = layer - 1:-1:0
                                                index_1 = lambda_offset(i_layer+1);
                                                index_2 = lambda_offset(i_layer+2);
                                                for beta = 0:index_1 - 1
                                                    P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                        sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), ...
                                                        abs(P(beta + index_1 + index_2)));
                                                    %                 operation_count_P=operation_count_P+1;

                                                end
                                            end
                                    end
                                    %-------------------------------------------------------
                                    %C = update_C(pac_params, i, C, u_esti(i+1));
                                    %-------------------------------------------------------
                                    phi_mod_2 = mod(i, 2);
                                    C(1, 1+phi_mod_2) = u_esti(i+1);

                                    if (phi_mod_2 == 1) && (i ~= N - 1)
                                        layer = bit_layer_vec(i + 1);
                                        for i_layer = 0:layer - 1
                                            index_1 = lambda_offset(i_layer+1);
                                            index_2 = lambda_offset(i_layer+2);
                                            for beta = index_1:2 * index_1 - 1
                                                C(beta + index_1, 2) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                                                C(beta + index_2, 2) = C(beta, 2);
                                                %             operation_count_C=operation_count_C+1;

                                            end
                                        end
                                        index_1 = lambda_offset(layer+1);
                                        index_2 = lambda_offset(layer+2);
                                        for beta = index_1:2 * index_1 - 1
                                            C(beta + index_1, 1) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                                            C(beta + index_2, 1) = C(beta, 2);
                                            %         operation_count_C=operation_count_C+1;
                                        end
                                    end
                                    %-------------------------------------------------------

                                end
                            else
                                %P = update_P(pac_params, i_start, P, C);
                                %-------------------------------------------------------
                                layer = llr_layer_vec(i_start + 1);
                                switch i_start %Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits.
                                    case 0
                                        index_1 = lambda_offset(n);
                                        for beta = 0:index_1 - 1
                                            P(beta + index_1) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                                            %             operation_count_P=operation_count_P+1;
                                        end
                                        for i_layer = n - 2:-1:0
                                            index_1 = lambda_offset(i_layer+1);
                                            index_2 = lambda_offset(i_layer+2);
                                            for beta = 0:index_1 - 1
                                                P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                %                 operation_count_P=operation_count_P+1;
                                            end
                                        end
                                    case N / 2
                                        index_1 = lambda_offset(n);
                                        for beta = 0:index_1 - 1
                                            x_tmp = C(beta + index_1,1);
                                            P(beta + index_1) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                                            %             operation_count_P=operation_count_P+1;

                                        end
                                        for i_layer = n - 2:-1:0
                                            index_1 = lambda_offset(i_layer+1);
                                            index_2 = lambda_offset(i_layer+2);
                                            for beta = 0:index_1 - 1
                                                P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
                                                %                 operation_count_P=operation_count_P+1;

                                            end
                                        end
                                    otherwise
                                        index_1 = lambda_offset(layer+1);
                                        index_2 = lambda_offset(layer+2);
                                        for beta = 0:index_1 - 1
                                            P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * P(beta + index_2) + ...
                                                P(beta + index_1 + index_2);
                                            %             operation_count_P=operation_count_P+1;

                                        end
                                        for i_layer = layer - 1:-1:0
                                            index_1 = lambda_offset(i_layer+1);
                                            index_2 = lambda_offset(i_layer+2);
                                            for beta = 0:index_1 - 1
                                                P(beta + index_1) = sign(P(beta + index_2)) * ...
                                                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), ...
                                                    abs(P(beta + index_1 + index_2)));
                                                %                 operation_count_P=operation_count_P+1;

                                            end
                                        end
                                end
                                %-------------------------------------------------------
                            end
                            %-----------------------------------------------------------

                            if delta(psi_temp+1) == 0
                                visited_before = true;
                                break;
                            elseif psi_temp == 0
                                Threshold = Threshold - Delta;
                                visited_before = false;
                                break;

                            end
                        else
                            Threshold = Threshold - Delta;
                            visited_before = false;
                            break;
                        end
                    end
                    psi = psi_temp;

                    %------------------------------------------------------------------------------
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
    d_esti(:,batch_index) = v_esti(info_bits_indices);
end
end
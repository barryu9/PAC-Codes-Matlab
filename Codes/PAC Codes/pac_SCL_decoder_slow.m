function d_esti = pac_SCL_decoder_slow(pac_params,rp, llr, L)
%LLR-based SCL deocoder, a single function, no other sub-functions.
%Frequently calling sub-functions will derease the efficiency of MATLAB
%codes.
%const

N = pac_params.N;
k = pac_params.k;
gen = pac_params.gen;

frozen_bits_mask = rp.frozen_bits_mask;
H_crc = pac_params.H_crc;
crc_length = pac_params.crc_length;
%memory declared
%If you can understand lazy copy and you just start learning polar codes
%for just fews days, you are very clever,
P = zeros(N - 1, L); %Channel llr is public-used, so N - 1 is enough.
C = zeros(N - 1, 2 * L); %I do not esitimate (x1, x2, ... , xN), so N - 1 is enough.
d = zeros(k, L); %unfrozen bits that polar codes carry, including crc bits.
PM = zeros(L, 1); %Path metrics
activepath = zeros(L, 1); %Indicate if the path is active. '1'→active; '0' otherwise.
cnt_u = 1; %information bit counter

%initialize
activepath(1) = 1;

curr_state = zeros(pac_params.conv_depth-1, L);
curr_state_temp = zeros(pac_params.conv_depth-1, L);
u_left = zeros(1, L);
u_right = zeros(1, L);

%decoding starts
%default: in the case of path clone, the origianl path always corresponds to bit 0, while the new path bit 1.
for phi = 0:N - 1
    for l_index = 1:L
        if activepath(l_index) == 0
            continue;
        end
        P(:, l_index) = update_P(pac_params, llr, phi, P(:, l_index), C(:, 2*l_index-1:2*l_index));
    end
    if frozen_bits_mask(phi + 1) == 0 %if now we decode an unfrozen bit
        PM_pair = realmax * ones(2, L);
        for l_index = 1:L
            if activepath(l_index) == 0
                continue;
            end
            curr_state_temp(:, l_index) = curr_state(:, l_index);
            [u_left(l_index), curr_state(:, l_index)] = conv_1bit_trans(0, curr_state(:, l_index), gen);
            [u_right(l_index), curr_state_temp(:, l_index)] = conv_1bit_trans(1, curr_state_temp(:, l_index), gen);
            PM_pair(1, l_index) = calc_PM(PM(l_index), P(1, l_index), u_left(l_index));
            PM_pair(2, l_index) = calc_PM(PM(l_index), P(1, l_index), u_right(l_index));
        end
        middle = min(2*sum(activepath), L);
        PM_sort = sort(PM_pair(:));
        PM_cv = PM_sort(middle);
        compare = PM_pair <= PM_cv;
        kill_index = zeros(L, 1); %to record the index of the path that is killed
        kill_cnt = 0; %the total number of killed path
        %the above two variables consist of a stack
        for i = 1:L
            if (compare(1, i) == 0) && (compare(2, i) == 0) %which indicates that this path should be killed
                activepath(i) = 0;
                kill_cnt = kill_cnt + 1; %push stack
                kill_index(kill_cnt) = i;
            end
        end
        for l_index = 1:L
            if activepath(l_index) == 0
                continue;
            end
            path_state = compare(1, l_index) * 2 + compare(2, l_index);
            switch path_state %path_state can equal to 0, but in this case we do no operation.
                case 1 % PM of the second row is lower
                    d(cnt_u, l_index) = 1;
                    C(:, 2*l_index-1:2*l_index) = update_C(pac_params, phi, C(:, 2*l_index-1:2*l_index), u_right(l_index));
                    PM(l_index) = PM_pair(2, l_index);
                    curr_state(:, l_index) = curr_state_temp(:, l_index);
                case 2 % PM of the first row is lower
                    d(cnt_u, l_index) = 0;
                    C(:, 2*l_index-1:2*l_index) = update_C(pac_params, phi, C(:, 2*l_index-1:2*l_index), u_left(l_index));
                    PM(l_index) = PM_pair(1, l_index);
                case 3 %
                    index = kill_index(kill_cnt);
                    kill_cnt = kill_cnt - 1; %pop stack
                    activepath(index) = 1;
                    %lazy copy
                    C(:, 2*index-1:2*index) = C(:, 2*l_index-1:2*l_index);
                    P(:, index) = P(:, l_index);
                    d(:, index) = d(:, l_index);
                    curr_state(:, index) = curr_state_temp(:, l_index);
                    d(cnt_u, l_index) = 0;
                    d(cnt_u, index) = 1;
                    C(:, 2*l_index-1:2*l_index) = update_C(pac_params, phi, C(:, 2*l_index-1:2*l_index), u_left(l_index));
                    C(:, 2*index-1:2*index) = update_C(pac_params, phi, C(:, 2*index-1:2*index), u_right(l_index));
                    PM(l_index) = PM_pair(1, l_index);
                    PM(index) = PM_pair(2, l_index);
            end
        end
        cnt_u = cnt_u + 1;
    else %frozen bit operation
        for l_index = 1:L
            if activepath(l_index) == 0
                continue;
            end
            [u_temp, curr_state(:, l_index)] = conv_1bit_trans(0, curr_state(:, l_index), gen);
            PM(l_index) = calc_PM(PM(l_index), P(1, l_index), u_temp);
            C(:, 2*l_index-1:2*l_index) = update_C(pac_params, phi, C(:, 2*l_index-1:2*l_index), u_temp);
        end
    end
end
%path selection.
activepath = logical(activepath);
PM_active = PM(activepath);
d_active = d(:, activepath);
[~, path_ordered] = sort(PM_active);
if (crc_length > 0)
    for l_index = 1:length(PM_active)
        path_num = path_ordered(l_index);
        info_with_crc = d(:, path_num);
        err = sum(mod(H_crc*info_with_crc, 2));
        if err == 0
            d_esti = d_active(1:end-crc_length, path_num);
            break;
        else
            if l_index == length(PM_active)
                d_esti = d_active(1:end-crc_length, path_ordered(1));
            end
        end
    end
else
    d_esti = d_active(1:end-crc_length, path_ordered(1));
end
end



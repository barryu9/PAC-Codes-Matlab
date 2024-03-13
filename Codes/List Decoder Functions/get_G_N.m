function G_N = get_G_N(F_N,N)
G_N = zeros(N, N);
G_N(1 : 2, 1 : 2) = F_N;
for i = 2 : log2(N)
    G_N(1 : 2^i, 1 : 2^i) = kron(G_N(1 : 2^(i - 1), 1 : 2^(i - 1)), F_N);
end

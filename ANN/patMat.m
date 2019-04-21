function ANN = patMat(i_A, L_A, n_A, f_A, i_B, L_B, n_B, f_B, radius)

[row, col] = size(L_A);
% find adjance patch and label it with 1
adj_A = find_adj_m(L_A);
adj_B = find_adj_m(L_B);

% list of barycenter of each patch
sA = regionprops(L_A, 'centroid');
c_A= cat(1, sA.Centroid);
sB = regionprops(L_B, 'centroid');
c_B = cat(1, sB.Centroid);

% find patch match within radius equal to 1
spm_A = gen_spm(c_A, radius);
spm_B = gen_spm(c_B, radius);

% patch match Inizializes as ANN(i) = i
ANN = (1:1:n_A)';

for i = 1:n_A-1

    for j = 1:4
        ANN = patMat_iter(i, j, ANN, adj_A, c_A, spm_A, n_A, f_A, ...
            adj_B, c_B, spm_B, n_B, f_B, radius, row, col);
    end
    %fprintf('%d ', i)    
end

disp('%n Done! ')

end
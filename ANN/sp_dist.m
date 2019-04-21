function Output = sp_dist(Ai, Bj, c_Ai, c_Bj,c_A, c_B, radius, nA, ...
    F_A, F_B, row, col)
% find super patch distance base on give radius
num_ij = 0;
den_ij = 0;
    
for k = Ai
    nc_Ak = c_A(k, :);
    nF_Ak = F_A(k);
    for l = Bj
        nc_Bl = c_B(l, :);
        nF_Bl = F_B(l);
           
        weight_ij = sp_weight(c_Ai, nc_Ak, c_Bj, nc_Bl, row, col, nA, radius);
           
        num_ij = num_ij + (weight_ij * norm(nF_Bl - nF_Ak));
        den_ij = den_ij + weight_ij;
    end
end

Output = num_ij / den_ij;
end


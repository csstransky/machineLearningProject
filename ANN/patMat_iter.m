function Output = SPM_iter(i, j, ANN, adj_A, c_A, spm_A, n_A, f_A, ...
    adj_B, c_B, spm_B, n_B, f_B, radius, row, col)

pooling = 0;
thresholdAngle = inf;
Ai = find(spm_A(i, :) == 1);
ANNi = find(spm_B(ANN(i), :) == 1);
    
% distance between Ai and its current ANN
center_Ai = c_A(i, :);
center_ANNi = c_B(ANN(i), :);
    
% sp_dist(i, ANN(i)) = min_dist;
min_dist = sp_dist(Ai, ANNi, center_Ai, center_ANNi, c_A, c_B, radius, ...
    n_A, f_A, f_B, row, col);
    
find_adj_Ai = find(adj_A(i, :) == 1);
find_adj_Ai = find_adj_Ai(find_adj_Ai ~= i);
    
% Dual interation
if (mod(j, 2) == 0)
    for i_ = find_adj_Ai
        anglei_i = sp_angle(c_A(i, :), c_A(i_, :));
        if (0.25*pi <= anglei_i && anglei_i < 1.25*pi)
            find_adj_Bi = find(adj_B(i_, :) == 1);
            find_adj_Bi = find_adj_Bi(find_adj_Bi ~= i_);
                
            for k = find_adj_Bi
                % add overflower checker
                if k>length(c_B)
                    break;
                end
                anglei_k = sp_angle(c_B(i_, :), c_B(k, :)); 
                absAngle = abs(anglei_k - anglei_i);
                if (absAngle < thresholdAngle)
                    pooling = k;
                    thresholdAngle = absAngle;
                end
            end
        end
            
        % varify the distance
        try
            superpatch_Bk = find(spm_B(pooling, :) == 1);
            center_Bk = c_B(pooling, :);
            comp_dist = sp_dist(Ai, superpatch_Bk, center_Ai, center_Bk, ...
                c_A, c_B, radius, n_A, f_A, f_B, row, col);
            % save the minimun distance super patch to ANN
            if (comp_dist < min_dist)
                ANN(i) = pooling;
                min_dist = comp_dist;
            end
        catch
        end
    end
        
% Single interation
else
    for i_ = find_adj_Ai
        anglei_i = sp_angle(c_A(i, :), c_A(i_, :));
        if (anglei_i < 0.25*pi || (1.25*pi <= anglei_i && anglei_i <= 2*pi))
            find_adj_Bi = find(adj_B(i_, :) == 1);
            find_adj_Bi = find_adj_Bi(find_adj_Bi ~= i_);

            for k = find_adj_Bi
                % add overflower checker
                if k>length(c_B)
                    break;
                end    
                anglei_k = sp_angle(c_B(i_, :), c_B(k, :)); 
                absAngle = abs(anglei_k - anglei_i);
                if (absAngle < thresholdAngle)
                    pooling = k;
                    thresholdAngle = absAngle;
                end
            end
        end
                
        % varify the distance
        try
            superpatch_Bk = find(spm_B(pooling, :) == 1);
            center_Bk = c_B(pooling, :);
            comp_dist = sp_dist(Ai, superpatch_Bk, center_Ai, center_Bk, ...
                c_A, c_B, radius, n_A, L_A, f_A, f_B, row, col);
            % save the minimun distance super patch to ANN
            if (comp_dist < min_dist)
                ANN(i) = pooling;
                min_dist = comp_dist;
            end
        catch
        end
    end
end
    
% find the number of superpixels to test from image B
testset = randperm(n_B, 8); 
for j = testset
   if j ~= ANN(i)
      superpatch_Bj = find(spm_B(j,:) == 1);
      comp_dist =  sp_dist(Ai, superpatch_Bj, center_Ai, c_B(j, :), ...
          c_A, c_B, radius, n_A, f_A, f_B, row, col);
      dist_m(i, j) = comp_dist;
      if (comp_dist < min_dist)
          ANN(i) = j;
          min_dist = comp_dist;
      end
   end
end
    
Output = ANN;

end


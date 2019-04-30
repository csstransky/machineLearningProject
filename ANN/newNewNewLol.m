close all; clear; clc;

% for reproducibility
rng('default');                             

cd('../trainImages');
train_img_files = dir('*.tif');
cd('../trainResults');
train_gt_files = dir('*.tif');
cd('../ANN')

% Cropped AerialImageDataset contains 18,001 images. Randomly choose N
N = 20;
radius = 20;
train_image = cell(N, 1);
ground_true = cell(N, 1);
X = randperm(length(train_img_files), N);
train_img_files = train_img_files(X);
train_gt_files = train_gt_files(X);
n = zeros(N, 1);


% patch-based approximate nearest neighbor (ANN) search methods are used to
% find correspondences. 
ANN = cell(N, N);

% Requested 5000x5000x100 (18.6GB) array exceeds maximum array size preference. 
% Creation of arrays greater than this limit may take a long time and cause 
% MATLAB to become unresponsive. 
% Data = zeros(5000, 5000, N);
Data = zeros(500, 500, N);
tic;

%% ANN initialization step.
disp('Step 1: ANN initialization step. ')
for i = 1:N
   img = imread(strcat('C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project\trainImages\', train_img_files(i).name));
   gt = imread(strcat('../trainResults/', train_gt_files(i).name));
   [l, h, w] = size(img);
   [L, nLables] = superpixels(img, floor(sqrt(l*h)));
   
   % ANN initialization step, For each superpixel Ai in A, we assign a 
   % random superpixel B(i) in B.
   L = patMat_init(L); 
   
   % image too large for 5000x5000
   Data(:, :, i) = L;
   n(i) = nLables;
   train_image{i} = img;
   ground_true{i} = gt;
   fprintf('%d ', i)
end
disp('%n Done! ')


%% ANN propagation step. 
disp('Step 2: ANN propagation step.')
labels = cell(N, 1);

% assign features to a variable
features = cell(N, 1); 
for i = 1:N
    labels{i} = label_img(ground_true{i}, Data(:, :, i), n(i));
    for j = 1:n(i)
        features{i} = [features{i} ; m_tri(j, Data(:, :, i), ...
            double(train_image{i}))];
    end
    fprintf('%d ', i)
end

disp('%n Done! ')          

%% ANN random search step. 
for i = 1:N
    i_A = train_image{i};
    L_A = Data(:, :, i);
    n_A = n(i);
    f_A = features{i};
    for j = (i+1):N
        fprintf('\nperformaning ANN on images %d and %d:\n', i, j)
        i_B = train_image{j};
        L_B = Data(:, :, j);
        n_B = n(j);
        f_B = features{j};
        
        % rearrange if the the image don't have enought superpixel
        if (n_B < n_A)
            temp=i_A; i_A=i_B; i_B=temp;
            temp=L_A; L_A=L_B; L_B=temp;
            temp=n_A; n_A=n_B; n_B=temp;
            temp=f_A; f_A=f_B; f_B=temp;
        end
        
        % Assign pathc match result to approximate nearest neighor matrix
        ANN{i, j} = patMat(i_A, L_A, n_A, f_A, i_B, L_B, n_B, f_B, radius); 
    end
end

% show total time spent on running step 1 to 3
fprintf('Total time spent: %f\n', toc)

%% Evaluate solution  

% Save ANN results
save('output','ANN','Data','features','ground_true','train_image',...
    'labels','n');

Loss = 0;
for l = 1:length(n)
    IMG_ = train_image{l};
    GT_ = logical(ground_true{l});
    
    threshold = imbinarize(rgb2gray(IMG_), graythresh(IMG_));
    Loss = Loss+mean(mean(abs(GT_-threshold)));
end

Loss = Loss/N;
fprintf('Loss of accuarcy during binary transformation = %f\n', Loss)

% find correctness
n_correct = ones(length(n)) .* Inf;
for i = 1:length(n)
    labels_A = labels{i};
    n_A = n(i);
    for j = (i+1):length(n)

        ANN_AB = ANN{i, j};
        n_B = n(j);
        try
            labels_B = labels{j};

            labels_B_est = ones(n_B, 1) .* inf;

            for k = 1:length(labels_A)
                labels_B_est(ANN_AB(k)) = labels_A(k);
            end
            
            l_compare = [labels_B, labels_B_est];
            l_compare = l_compare(labels_B_est ~= Inf, :);
            count_l = l_compare(:, 1) == l_compare(:, 2);
            n_correct(i, j) = length(count_l(count_l == 1))/length(count_l);
        catch      
        end
    end
end

final_result = mean(mean(n_correct(n_correct ~= Inf)));

fprintf('Percentage correctness = %f\n', final_result*100)

function Output = sp_weight(c_Ai, c_Ai_, c_Bj, c_Bj_, m, n, nA, radius)

    ws_Ai_ = exp(-((c_Ai(2)-c_Ai_(2))^2+(c_Ai(1)-c_Ai_(1))^2)/(2*radius^2));
    ws_Bj_ = exp(-((c_Bj(2)-c_Bj_(2))^2+(c_Bj(1)-c_Bj_(1))^2)/(2*radius^2));

    Output = exp(-(c_Bj_ - c_Ai_+ c_Ai - c_Bj)*(c_Bj_-c_Ai_+c_Ai-c_Bj)'/(0.5 ...
        *sqrt(m*n/nA))) * ws_Ai_ * ws_Bj_;
end

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

function Output = sp_angle(c_A, c_B)
    % find super patch angle base on given centers barycenters
    x = c_B(1) - c_A(1);
    y = c_B(2) - c_A(2);

    angle = atan2(y, x);

    if (c_B(2) < c_A(2))
        Output = -angle;
    else
        Output = 2*pi - angle;
    end
end

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

function Output = patMat_init(X)
    Output = zeros(size(X));
    Y = randperm(max(max(X)));

    for i = 1:length(Y)
        Output(find(X == i)) = Y(i);
    end
end

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

function Output = m_tri(idx, Data, image)
    % Calculate the mean of the tristimulus values of Data..

    red_value = (Data == idx) .* image(:, :, 1);
    green_value = (Data == idx) .* image(:, :, 2);
    blue_value = (Data == idx) .* image(:, :, 3);

    RED = mean(mean(red_value(find(red_value ~= 0))));
    GREEN = mean(mean(green_value(find(green_value ~= 0))));
    BLUE = mean(mean(blue_value(find(blue_value ~= 0))));

    Output = [RED, GREEN, BLUE];
end

function Output = gen_spm(center, radius)
    % find the super patches fom the given radius
    Output = zeros(length(center));

    for i = 1:length(center(:, 1))
        temp1 = center(i, :);

        for j = 1:length(center(:, 1))
            temp2 = center(j, :);
            if (temp2(1)-temp1(1))^2+(temp2(2)-temp1(2))^2<=radius^2
                Output(i, j) = 1;
            end
        end
    end
end

function Output = find_adj_m(X)
    [row, col]=size(X);
    Output=zeros(row, col);

    for i=1:row
        for j=1:col
            if abs(i-j)==1
                Output(i,j)=1;
            end
        end
    end
end

function Output = label_img(gt, Data, numIdx)
    % label Data using the barycenter of superpatch of gt

    Y = regionprops(Data, 'centroid');
    barycenter = round(cat(1, Y.Centroid));
    Output = zeros(numIdx, 1);

    for i = 1:numIdx
       Output(i) = gt(barycenter(i, 2), barycenter(i, 1)); 
    end
end


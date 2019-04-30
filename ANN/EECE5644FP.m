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



close all; clear all; clc

% for reproducibility
rng('default');

cd('C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project\trainImages\');
train_img_files = dir('*.tif');
cd('C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project\trainResults\');
train_gt_files = dir('*.tif');
cd('/..')

% Cropped AerialImageDataset contains 18,001 images. Randomly choose N
N = 20;
radius = 20;
train_image = cell(N, 1);
ground_true = cell(N, 1);
X = randperm(length(train_img_files), N);
train_img_files = train_img_files(X);
train_gt_files = train_gt_files(X);
n = zeros(N, 1);
K = 8;

tic;
accuracy = [];
for i = 1:N
   img = im2double(imread(strcat('C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project\trainImages\', train_img_files(i).name)));
   gt = im2double(imread(strcat('C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project\trainResults\', train_gt_files(i).name)));
   new_I = reshape(img,size(img,1)*size(img,2),3);

    c_center = new_I( ceil(rand(K,1)*size(new_I,1)) ,:);
    % Distances and Labels
    DAL = zeros(size(new_I,1),K+2);
    % K-means Iteration
    KMI = 8;
    for n = 1:KMI
       for i = 1:size(new_I,1)
           for j = 1:K
               DAL(i,j) = norm(new_I(i,:) - c_center(j,:));
           end
           % 1:K are Distance from Cluster Centers 1:K
           [Distance, labels] = min(DAL(i,1:K));
           % K+1 is Cluster Label
           DAL(i,K+1) = labels;
           % K+2 is Minimum Distance
           DAL(i,K+2) = Distance;
       end
       for i = 1:K
           % Cluster K Points
           A = (DAL(:,K+1) == i);
           % New Cluster Centers
           c_center(i,:) = mean(new_I(A,:));
           % If CENTS(i,:) Is Nan Then Replace It With Random Point
           if sum(isnan(c_center(:))) ~= 0
               % Find Nan Centers
               NC = find(isnan(c_center(:,1)) == 1);
               for Ind = 1:size(NC,1)
               c_center(NC(Ind),:) = new_I(randi(size(new_I,1)),:);
               end
           end
       end
    end

    X = zeros(size(new_I));

    for i = 1:K
       idx = find(DAL(:,K+1) == i);
       X(idx,:) = repmat(c_center(i,:),size(idx,1),1);
    end

    %seg_I = reshape(X,size(I,1),size(I,2),3);
    seg_I = round(rgb2gray(reshape(X,size(img,1),size(img,2),3)));

    match=0;
    [row, col, hig] = size(img);
    for i=1:row
       for j=1:col
           if gt(i,j) == seg_I(i,j)
               match=match+1;
           end
       end
    end
    match=match/(row*col);
    accuracy = [accuracy;match];
    fprintf('correctness = %f\n', match*100)
end
toc
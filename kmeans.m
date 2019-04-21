close all; clear; clc;% for reproducibility
rng('default');% Read an image into the workspace.
cd('trainImages/');
I = im2double(imread('austin1_cropped_28.tif'));
cd('../trainResults/');
gt = im2double(imread('austin1_cropped_28.tif'));
cd('..')

%imshow(I)

%title(‘Original Image’)
new_I = reshape(I,size(I,1)*size(I,2),3);
%Segment the image into three regions using k-means clustering.
%[L,Centers] = imsegkmeans(I,3);
%B = labeloverlay(I,L);
%imshow(B)

%title(‘Labeled Image’)% Cluster Numbers
K = 8;

% Cluster Centers
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
seg_I = round(rgb2gray(reshape(X,size(I,1),size(I,2),3)));

%% Test accuracy
match=0;
[row, col, hig] = size(I);
for i=1:row
   for j=1:col
       if seg_I(i,j)==gt(i,j)
           match=match+1;
       end
   end
end

match=match/(row*col);
fprintf('correctness = %f\n', match*100)
imshow(seg_I)
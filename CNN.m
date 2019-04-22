
training = false; % SET THIS TO "true" TO TRAIN THE NETWORK
if ~training
    load('cnn_parameters.mat');
end
    
%localFolder = 'C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project';
localFolder = '/home/bazooka/Documents/EECE5644_Machine_Learn_Pattern_Recogntn/project';
imageFolder = fullfile(localFolder, 'trainImagesSelect');
resultsFolder = fullfile(localFolder, 'trainResultsSelect');
imds = imageDatastore(imageFolder, 'LabelSource', 'foldernames', 'IncludeSubfolders',true);

if training
    buildingRGB = 255;
    terrainRGB = 0;
    classNames = ["buildings", "terrain"];
    pixelLabelID = [buildingRGB terrainRGB];
    pxds = pixelLabelDatastore(resultsFolder,classNames,pixelLabelID);

    %{
    C = readimage(pxds,3);
    C(5,200)

    I = readimage(imds,3);
    I = histeq(I);
    imshow(I);

    B = labeloverlay(I,C);
    figure
    imshow(B)
    %}

    I = read(imds);

    numFilters = 8;
    filterSize = 3;
    numClasses = 2;
    layers = [
        imageInputLayer([500 500 3])
        convolution2dLayer(filterSize,numFilters,'Padding',1)
        reluLayer()
        maxPooling2dLayer(2,'Stride',2)
        convolution2dLayer(filterSize,numFilters,'Padding',1)
        reluLayer()
        transposedConv2dLayer(4,numFilters,'Stride',2,'Cropping',1);
        convolution2dLayer(1,numClasses);
        softmaxLayer()
        pixelClassificationLayer()
        ];

    opts = trainingOptions('sgdm', ...
        'InitialLearnRate',0.001, ...
        'MaxEpochs',2, ...
        'MiniBatchSize',8);

    trainingData = pixelLabelImageDatastore(imds,pxds);


    tbl = countEachLabel(trainingData);
    totalNumberOfPixels = sum(tbl.PixelCount);
    frequency = tbl.PixelCount / totalNumberOfPixels;
    classWeights = 1./frequency;
    layers(end) = pixelClassificationLayer('Classes',tbl.Name,'ClassWeights',classWeights);
    net = trainNetwork(trainingData,layers,opts);
end

testImage = readimage(imds,2);
imshow(testImage);
C = semanticseg(testImage,net);
B = labeloverlay(testImage,C);
imshow(B)


% Images are 500x500

% 20 images took an 1.5 hours for neariest neighbor
% time; performance; error
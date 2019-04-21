% Location of folder on home computer
localFolder = 'C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project';
imageFolder = fullfile(localFolder, 'trainImagesSelect');
resultsFolder = fullfile(localFolder, 'trainResultsSelect');
imds = imageDatastore(imageFolder, 'LabelSource', 'foldernames', 'IncludeSubfolders',true);

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

%{
% Create an image input layer
inputSize = [32 32 3];
imgLayer = imageInputLayer(inputSize);

% Create Downsampling Network
filterSize = 3;
numFilters = 32;
conv = convolution2dLayer(filterSize,numFilters,'Padding',1);
relu = reluLayer();

poolSize = 2;
maxPoolDownsample2x = maxPooling2dLayer(poolSize,'Stride',2);
downsamplingLayers = [
    conv
    relu
    maxPoolDownsample2x
    conv
    relu
    maxPoolDownsample2x
    ];

% Create upsampling Network
filterSize = 4;
transposedConvUpsample2x = transposedConv2dLayer(4,numFilters,'Stride',2,'Cropping',1);
upsamplingLayers = [
    transposedConvUpsample2x
    relu
    transposedConvUpsample2x
    relu
    ];

% Create a pixel classification layer
numClasses = 3;
conv1x1 = convolution2dLayer(1,numClasses);
finalLayers = [
    conv1x1
    softmaxLayer()
    pixelClassificationLayer()
    ];

% Complete network layer stack
net = [
    imgLayer    
    downsamplingLayers
    upsamplingLayers
    finalLayers
    ]
%}

I = read(imds);
C = read(pxds);


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


%{
net = resnet18;
%net = trainNetwork(trainingData,layers,opts);
testImage = readimage(imds,2);
imshow(testImage);
C = semanticseg(testImage,net);
B = labeloverlay(testImage,C);
imshow(B);
%}

tbl = countEachLabel(trainingData);
totalNumberOfPixels = sum(tbl.PixelCount);
frequency = tbl.PixelCount / totalNumberOfPixels;
classWeights = 1./frequency;
layers(end) = pixelClassificationLayer('Classes',tbl.Name,'ClassWeights',classWeights);
net = trainNetwork(trainingData,layers,opts);
testImage = readimage(imds,2);
imshow(testImage);
C = semanticseg(testImage,net);
B = labeloverlay(testImage,C);
imshow(B)


% Images are 500x500

% 20 images took an 1.5 hours for neariest neighbor
% time; performance; error
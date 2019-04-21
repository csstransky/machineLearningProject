pretrainedURL = 'https://www.mathworks.com/supportfiles/vision/data/deeplabv3plusResnet18CamVid.mat';
pretrainedFolder = fullfile(tempdir,'pretrainedNetwork');
pretrainedNetwork = fullfile(pretrainedFolder,'deeplabv3plusResnet18CamVid.mat'); 
if ~exist(pretrainedFolder,'dir')
    mkdir(pretrainedFolder);
    disp('Downloading pretrained network (58 MB)...');
    websave(pretrainedNetwork,pretrainedURL);
end

% Location of folder on home computer
localFolder = 'C:\Users\Bazooka\Desktop\EECE5644_Machine_Learn_Pattern_Recogntn\project';
imageFolder = fullfile(localFolder, 'trainImages');
resultsFolder = fullfile(localFolder, 'trainResults');
imds = imageDatastore(imageFolder, 'LabelSource', 'foldernames', 'IncludeSubfolders',true);

buildingRGB = 255;
terrainRGB = 0;
classes = ["buildings", "terrain"];
pixelLabelID = [buildingRGB terrainRGB];
pxds = pixelLabelDatastore(resultsFolder,classes,pixelLabelID);


C = readimage(pxds,3);
C(5,200)

I = readimage(imds,3);
I = histeq(I);
imshow(I);

B = labeloverlay(I,C);
figure
imshow(B)
[imdsTrain, imdsVal, imdsTest, pxdsTrain, pxdsVal, pxdsTest] = partitionCamVidData(imds,pxds);
numTrainingImages = numel(imdsTrain.Files);
numValImages = numel(imdsVal.Files);
numTestingImages = numel(imdsTest.Files);

% Specify the network image size. This is typically the same as the traing image sizes.
imageSize = [500 500 3];

% Specify the number of classes.
numClasses = numel(classes);

% Create DeepLab v3+.
lgraph = helperDeeplabv3PlusResnet18(imageSize, numClasses);
tbl = countEachLabel(pxds);
imageFreq = tbl.PixelCount ./ tbl.ImagePixelCount;
classWeights = median(imageFreq) ./ imageFreq;
pxLayer = pixelClassificationLayer('Name','labels','Classes',tbl.Name,'ClassWeights',classWeights);
lgraph = replaceLayer(lgraph,"classification",pxLayer);

% Define validation data.
pximdsVal = pixelLabelImageDatastore(imdsVal,pxdsVal);

% Define training options. 
options = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',10,...
    'LearnRateDropFactor',0.3,...
    'Momentum',0.9, ...
    'InitialLearnRate',1e-3, ...
    'L2Regularization',0.005, ...
    'ValidationData',pximdsVal,...
    'MaxEpochs',30, ...  
    'MiniBatchSize',8, ...
    'Shuffle','every-epoch', ...
    'CheckpointPath', tempdir, ...
    'VerboseFrequency',2,...
    'Plots','training-progress',...
    'ValidationPatience', 4); ...
    
augmenter = imageDataAugmenter('RandXReflection',true,...
    'RandXTranslation',[-10 10],'RandYTranslation',[-10 10]);
pximds = pixelLabelImageDatastore(imdsTrain,pxdsTrain, ...
    'DataAugmentation',augmenter);

doTraining = true;
if doTraining    
    [net, info] = trainNetwork(pximds,lgraph,options);
else
    %data = load(pretrainedNetwork); 
    %net = data.net;
end

I = readimage(imdsTest,10);
C = semanticseg(I,net);
B = labeloverlay(I,C);
imshow(B)

function buildingDetection(inputFileName, outputFileName)

[imgFile, lidarFile, demFile, LPDx, LPDy, useEntropy, useAdjustment, useRefinement, useReduction] = readInputFile(inputFileName);

out = 'reading data ...'
[A, R, bbox] = geotiffread(imgFile);
[Ad, Rd, bboxd] = geotiffread(demFile);
fp = fopen(lidarFile, 'r');
ALS = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
ALS = ALS';

out = 'generating masks ...'
[pMask sMask Rm] = generateMask(A, R, Ad, Rd, ALS, LPDx, LPDy);

out = 'extracting lines from the primary mask ...'
extractLines(pMask, Rm, 'maskLines.txt');

out = 'computing NDVI ...'
[ndvi ndvisig] = computeNDVI(A, 0);%second argument: 1 to show ndvi image; 0 otherwise

[height width bands] = size(A);
if bands == 3 %we need sigma of pseudo NDVI, so convert ndvisig to 8-bit greyscale image
    mn = min(min(ndvisig));
    mx = max(max(ndvisig));
    ndvisig = uint8((ndvisig-mn)*255/(mx-mn)+0.5);    
    ndviThresh = 48;
    ndvi = ndvisig;
else
    ndviThresh = 10;
end


if useEntropy  
    out = 'computing image entropy ...'
    % calculate image texture
    Ig = rgb2gray(A(:,:,1:3));
    E = entropyfilt(Ig);
    Eim = mat2gray(E);    
    eMask = im2bw(Eim, .8);
    %figure; imshow(eMask);
    useRefinement = 1;
else
    eMask = [];
    %useRefinement = 0;
end

if (useAdjustment || useRefinement)
    out = 'extracting lines and/or orientation information from the image ...'
    [orinetation edgeimage] = extractImageLines(A,R,useAdjustment,useRefinement,'imageLines.txt');
end

out = 'applying NDVI/entropy thresholds ...'
applyNDVIentropy(A, R, pMask, Rm, ndvi, ndviThresh, eMask, useEntropy, useAdjustment, 'maskLines.txt', 'ndviLines.txt');

out = 'extending lines ...'
extendLines(A, R, bbox, pMask, Rm, ndvi, ndviThresh, eMask, useEntropy, 'ndviLines.txt', 'extendedLines.txt');

out = 'forming initial candidates ...'
formCandidates(A, R, bbox, pMask, Rm, ndvi, ndviThresh, eMask, useEntropy, 'extendedLines.txt', 'initialCandidates.txt');

out = 'extending candidates ...'
extendCandidates(A, R, pMask, sMask, Rm, useReduction, 'initialCandidates.txt', 'finalCandidates.txt')

if useRefinement
    out = 'refining candidates ...'
    refineCandidates(A, R, orinetation, edgeimage, 'finalCandidates.txt', outputFileName);
end



function [imgFile, lidarFile, demFile, LPDx, LPDy, useEntropy, useAdjustment, useRefinement, useReduction] = readInputFile(inputFileName)
fp = fopen(inputFileName,'r');

%read image file name
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
imgFile = ins;

%read lidar file name
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
lidarFile = ins;

%read dem file name
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
demFile = ins;

%read LIDAR point-to-point distance in x direction (in metre)
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
LPDx = str2double(ins);

%read LIDAR point-to-point distance in y direction (in metre)
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
LPDy = str2double(ins);

%read useEntropy flag
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
useEntropy = str2double(ins);

%read useAdjustment flag
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
useAdjustment = str2double(ins);

%read useRefinement flag
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
useRefinement = str2double(ins);

%read useRefinement flag
ins = fgetl(fp);
while (size(ins,2) < 1 || strcmp(ins(1,1),'%'))
    ins = fgetl(fp);
end
useReduction = str2double(ins);
fclose(fp);

%% Look for Centers of Potential Clusters
clearvars; close all; clc;
I = rgb2gray(imread('E:\Download\Spring 2019\Machine Learning\EECE5644_ANN\images\austin21_cropped_37.tif'));
level = multithresh(I);
seg_I = imquantize(I,level);
seg_I = seg_I - 1;
SE = strel('octagon',6);
J = imerode(seg_I,SE);
J = J > 0;
stats = regionprops(J,I,{'Centroid','MajorAxisLength','MinorAxisLength'});

figure(2)
imshow(I)
hold on
for k = 1:numel(stats)
    plot(stats(k).Centroid(1),stats(k).Centroid(2),'r*')
end
hold off

% Some Thresholds
diagonal_threshold = 7;
radiusMargin1 = 5;
radiusMargin2 = 1;
horizontalAlignment1 = 5;
horizontalAlignment2 = 10;
verticalTolerance1 = 0.3;
verticalTolerance2 = 0.1;
theta90_threshold = 5;
%{
%% Edge Detection with Fuzzy Logic
sharpened_I = imsharpen(I,'radius',2,'Amount',1);
tempI = im2double(sharpened_I);
Gx = [-1 1];
Gy = Gx';
Ix = conv2(tempI,Gx,'same');
Iy = conv2(tempI,Gy,'same');
edgeFIS = mamfis('Name','edgeDetection');
edgeFIS = addInput(edgeFIS,[-1 1],'Name','Ix');
edgeFIS = addInput(edgeFIS,[-1 1],'Name','Iy');
sx = 0.1;
sy = 0.1;
edgeFIS = addMF(edgeFIS,'Ix','gaussmf',[sx 0],'Name','zero');
edgeFIS = addMF(edgeFIS,'Iy','gaussmf',[sy 0],'Name','zero');
edgeFIS = addOutput(edgeFIS,[0 1],'Name','Iout');
wa = 0.1;
wb = 1;
wc = 1;
ba = 0;
bb = 0;
bc = 0.7;
edgeFIS = addMF(edgeFIS,'Iout','trimf',[wa wb wc],'Name','white');
edgeFIS = addMF(edgeFIS,'Iout','trimf',[ba bb bc],'Name','black');
r1 = "If Ix is zero and Iy is zero then Iout is white";
r2 = "If Ix is not zero or Iy is not zero then Iout is black";
edgeFIS = addRule(edgeFIS,[r1 r2]);
Ieval = zeros(size(I));
for ii = 1:size(I,1)
    Ieval(ii,:) = evalfis(edgeFIS,[(Ix(ii,:));(Iy(ii,:))]');
end
%}

%% Edge Detection with Sobel
sharpened_I = imsharpen(I,'radius',2,'Amount',1);
[~,threshold] = edge(sharpened_I,'sobel');
fudgeFactor = 0.5;
bw = edge(sharpened_I,'sobel',threshold * fudgeFactor);

%% Mark Label
finalBW = poly2mask(0,0,500,500);
for i = 8
    diagonal = (sqrt(stats(i).MajorAxisLength^2 + stats(i).MinorAxisLength^2));
    if diagonal >= diagonal_threshold
        center = (stats(i).Centroid);
        outerR = diagonal;
        innerR = (stats(i).MinorAxisLength);
        outerR = outerR/2 + radiusMargin1;
        innerR = innerR/2 - radiusMargin2;
        %thisWindow = bw(center(2)-outerR:center(2)+outerR, center(1)-outerR:center(1)+outerR);
        thisWindow = imcrop(bw,[center(1)-outerR center(2)-outerR diagonal+2*radiusMargin1-1 diagonal+2*radiusMargin1-1]);
        [xgrid, ygrid] = meshgrid(1:size(thisWindow,2),1:size(thisWindow,1));
        ringMask = ((xgrid-size(thisWindow,2)/2).^2 + (ygrid-size(thisWindow,1)/2).^2) <= outerR.^2 & ((xgrid-size(thisWindow,2)/2).^2 + (ygrid-size(thisWindow,1)/2).^2) >= innerR^2;
        thisRing = thisWindow .* ringMask;
        [H,T,R] = hough(thisRing,'RhoResolution',1,'Theta',-90:89.5);
        
        % Split H and T in halves along the theta axis and find local maxima separately
        left_H = H(:,1:90);
        right_H = H(:,91:180);
        left_T = T(:,1:90);
        right_T = T(:,91:180);
        
        % find local maxima in hough transform
        left_P = houghpeaks(left_H,10,'threshold',ceil(0.4*max(left_H(:))));
        right_P = houghpeaks(right_H,10,'threshold',ceil(0.4*max(right_H(:))));
        realPeaks1 = [];
        maxIntensity = 0;
        for j = 1:size(right_P,1)
            for k = 1:size(right_P,1)
                lm_A = [right_T(right_P(j,2)) R(right_P(j,1))];
                lm_B = [right_T(right_P(k,2)) R(right_P(k,1))];
                result = check_vertical_alignment(lm_A, lm_B, size(thisWindow,1), horizontalAlignment1, verticalTolerance1);
                avgIntensity1 = (right_H(right_P(j,1),right_P(j,2)) + right_H(right_P(k,1),right_P(k,2))) / 2;
                if result && (avgIntensity1 > maxIntensity)
                    realPeaks1 = [lm_A;lm_B];
                    maxIntensity = avgIntensity1;
                end
            end
        end
        
        realPeaks2 = [];
        maxIntensity2 = 0;
        for m = 1:size(left_P,1)
            for n = 1:size(left_P,1)
                lm_C = [left_T(left_P(m,2)) R(left_P(m,1))];
                lm_D = [left_T(left_P(n,2)) R(left_P(n,1))];
                result = check_vertical_alignment2(lm_C, lm_D, size(thisWindow,1), horizontalAlignment2, verticalTolerance2);
                avgIntensity2 = (left_H(left_P(m,1),left_P(m,2)) + left_H(left_P(n,1),left_P(n,2))) / 2;
                if result && (avgIntensity2 > maxIntensity2)
                    realPeaks2 = [lm_C;lm_D];
                    maxIntensity2 = avgIntensity2;
                end
            end
        end
        if ~isempty(realPeaks1) && ~isempty(realPeaks2)
            if abs(abs((realPeaks1(1,1)+realPeaks1(2,1))/2 - (realPeaks2(1,1)+realPeaks2(2,1))/2) - 90) <= theta90_threshold
                coords = drawRectangleonImageAtAngle(I,center',abs(realPeaks1(1,2)-realPeaks1(2,2)),abs(realPeaks2(1,2)-realPeaks2(2,2)),-abs((realPeaks1(1,1)+realPeaks1(2,1))/2));
                finalBW = poly2mask(coords(1,:),coords(2,:),500,500) + finalBW;
            end
        end
    end
end
imshow(finalBW)

figure(2)
%imshow(imadjust(rescale(H)),'XData',T,'YData',R,'InitialMagnification','fit');
H = right_H;
T = right_T;
P = right_P;
imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
axis on, axis normal; hold on;
%P = houghpeaks(H(:,end/2:end),10,'threshold',ceil(0.4*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');
%plot(T(peaks(:,2)),R(peaks(:,1)),'s','color','white');

%%
figure(3)
%imshow(imadjust(rescale(H)),'XData',T,'YData',R,'InitialMagnification','fit');
H = left_H;
T = left_T;
P = left_P;
imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
axis on, axis normal; hold on;
%P = houghpeaks(H(:,end/2:end),10,'threshold',ceil(0.4*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');
%plot(T(peaks(:,2)),R(peaks(:,1)),'s','color','white');

function decision = check_vertical_alignment(pointA, pointB, windowSize, threshold, tolerance)

    horzDifference = abs(pointA(1)-pointB(1));
    horzCheck = horzDifference <= threshold;
    vertDistance = abs((pointA(2)+pointB(2))/2);
    if horzCheck
        thisTheta = abs((pointA(1)+pointB(1))/2);
        temp = (vertDistance*cosd(thisTheta) - windowSize/2) / tand(thisTheta);
        geoCheck = abs(temp - (windowSize/2 - vertDistance*sind(thisTheta))) / temp <= tolerance;
        if geoCheck && (vertDistance > windowSize/2) && (abs(pointA(2)-pointB(2)) > 7)
            decision = 1;
        else
            decision = 0;
        end
    else
        decision = 0;
    end
end

function decision = check_vertical_alignment2(pointA, pointB, windowSize, threshold, tolerance)

    horzDifference = abs(pointA(1)-pointB(1));
    horzCheck = horzDifference <= threshold;
    vertDistance = abs((pointA(2)+pointB(2))/2);
    if horzCheck
        thisTheta = abs((pointA(1)+pointB(1))/2);
        temp = (windowSize/2 - vertDistance/sind(thisTheta)) * tand(thisTheta);
        geoCheck = abs(temp - windowSize/2) / temp <= tolerance;
        if geoCheck && (vertDistance < windowSize/2) && (abs(pointA(2)-pointB(2)) > 7)
            decision = 1;
        else
            decision = 0;
        end
    else
        decision = 0;
    end
end
function [ret retInd H] = validate_candidate_building(A,R,orinetation, edgeimage,P1,P2,P3,P4,minBinHeight,minAggregateBinHeight)
%ret = 0;
binDist = 5;
hSize = 180/binDist;
%minBinHeight = 30;
%{
minBinHeight1 = 90; % 6m length; for farifield: minBinHeight = 40
minBinHeight2 = 60; % 6m length; for farifield: minBinHeight = 40
minStdDeviation = 12;
minPerOfMaxPeak = 0.15;
minMax2MeanRatio = 4.0;
%}
%find edge slopes
%G = find_ALS(A,R,grad,P1,P2,P3,P4);

P1i = obj2obs(P1,R);
P2i = obj2obs(P2,R);
P3i = obj2obs(P3,R);
P4i = obj2obs(P4,R);

H = zeros(1,hSize);
Hdummy = zeros(1,hSize);
Hindall = [];

OriEdge(:,:,1) = orinetation;
OriEdge(:,:,2) = edgeimage;
oe = findImagePoints(OriEdge, P1i, P2i, P3i, P4i);
cNumbers = [];
curNum = 1;

for i = 1:size(oe,1)
    curveNum = oe(i,2);
    if curveNum > 212
        here = 1;
    end
    if (curveNum > 0)
        found = false;
        for j = 1:size(cNumbers,1)
            if cNumbers(j,1) == curveNum
                found = true;
                break;
            end
        end
        if found % curve already added and a histogram was created
        else %form a new histogram
            cNumbers(curNum,1) = curveNum;
            curNum = curNum+1;
            Hindall = [Hindall; Hdummy];
            j = size(Hindall,1);
        end
        theta = 180*oe(i,1)/pi;        
        theta = theta + 90; % range: -90 to 90 is transferred to 0 to 180
        index = ceil(theta/binDist);
        if (index <= 0)
            index = 1;
        elseif (index > hSize)
            index = hSize;
        end
        H(1,index) = H(1,index)+1;
        Hindall(j,index) = Hindall(j,index) + 1;
    end
end

if size(Hindall,1) > 0
    indSize = sum(Hindall,2);
    [mx mxid] = max(indSize);
    if mx > minBinHeight
        Hind = Hindall(mxid,:);
    else
        Hind = Hdummy;
    end
else
    Hind = Hdummy;
end

ret = evaluateHistogram(H,binDist,minBinHeight,minAggregateBinHeight);
retInd = evaluateHistogram(Hind,binDist,minBinHeight,minAggregateBinHeight);

%{
%find histogram
H = findHistogram(G,binDist);
Hind = findIndividualHistogram(G,binDist);

%validate this building
%ret = evaluateHistogram(H,binDist,minBinHeight2,minStdDeviation,minPerOfMaxPeak);
ret = evaluateHistogram(H,binDist,minBinHeight,minAggregateBinHeight);


retInd = [];
%if (ret(1,2) >= minBinHeight1 || ret(1,4)>= minMax2MeanRatio)
%else
    for i = 1:size(Hind,1)
        if (sum(Hind(i,:)) >= minBinHeight)
            r = evaluateHistogram(Hind(i,:),binDist,minBinHeight,minAggregateBinHeight);
            retInd = [retInd; r];
        %else
        %    retInd = [retInd; [NaN NaN NaN NaN NaN NaN NaN]];
        end    
    end
    %ret = retInd;    
%end
if (size(retInd,1)>1)
    [mx id] = max(retInd(:,2));
    retInd = retInd(id,:);
end
%here = 1;
%}
%function ret = evaluateHistogram(H,binDist,minBinHeight,minStdDeviation,minPerOfMaxPeak)
function ret = evaluateHistogram(H,binDist,minBinHeight,minAggregateBinHeight)
%ret = [];
[gau W] = find_Gaussian(0.5);
Y = H';
L = size(Y,1);
Y1 = [ones(W,1)*2*Y(1)-Y(W+1:-1:2);Y;ones(W,1)*2*Y(L)-Y(L-1:-1:L-W)];
YY = conv(Y1,gau);
Ys = YY(2*W+1:L+2*W);
%minBuildingLength = minBinHeight;

extremum = find_extremum01(Ys);
extAngle = (extremum*binDist)';
Yext = Ys(extremum,1);
[Yexts ids] = sort(Yext,'descend');
mn = mean(Ys);

v = [Yexts - mn];
%    Val(i,1) = sum(v{i}(v{i} > 0))/sum(Ys);
%    Val1(i,1) = sum(v{i}(v{i} > 0))/sum(Yext);
    %Val2(i,1) = sum(v{i}(v{i} > 0)); %standard deviation
    %Val3(i,1) = sum(Y);
    Val4 = [mean(v(v > 0)) mn];
    %Val5(i,1) = max(v{i});
    if (size(extremum,1)>0)
        if (extremum(1,1) == 1 && extremum(1,end) == size(H,2) && H(1,1) >= minBinHeight && H(1,end) >= minBinHeight)
            Val6 = H(1,1)+H(1,end);
        else
            Val6 = max(Y); %maximum
        end
    else
        Val6 = max(Y); %maximum
    end
    perMax = Val6/sum(H);
    Val7 = [];   
    for j = 1:size(v,1)
        if (v(j,1)>0)
            id = ids(j,1);
            Val7 = [Val7 extAngle(id,1)];
        end
    end
    subdiv = [];
    for j=1:size(Val7,2)-1
        sub = abs(Val7(1,j+1:end) - Val7(1,j));
        subdiv = [subdiv sub/90];                
    end
    %Val8 = 0;
    Val8 = sum(subdiv == 1) + sum(subdiv == 2);% + sum(subdiv{i} == 3) + sum(subdiv{i} == 4);
    
    %if (Val8 > 0 || Val6 > minBinHeight || Val4(1,1) > minStdDeviation || perMax > minPerOfMaxPeak)
    %Confidence = [Val8 > 0 Val6 > minBinHeight Val4(1,1) > minStdDeviation perMax > minPerOfMaxPeak];
    %ret = sum(Confidence)/size(Confidence,2);    
    %ret = [Val8 > 0 Val6 > minBinHeight Val4(1,1) > minStdDeviation perMax > minPerOfMaxPeak];
    
    
    extremum = find_extremum01(Y);
    Yext = Y(extremum,1);
    ids = find(Yext >= minAggregateBinHeight);
    if (size(ids,1) > 1)
        sm = 0;
        for i=1:size(ids,1)-1
            st = extremum(1,ids(i,1))+1;
            en = extremum(1,ids(i+1,1))-1;
            if (st <= en)
                sm = sm + mean(Y(st:en));
            else %in the case of successive maxima in which case this is not a building
                sm = sm+minAggregateBinHeight;
            end
        end
        av = sm/(size(ids,1)-1);
    else
        av = -1;
    end
    %consider wrap around histogram
    ret = [Val8 Val6 mn Val6/mn Val4(1,1) perMax sum(H) av];
    %here = 1;
    %end

%{    
    
function H = findHistogram(G,binDist)
    H = zeros(1, 180/binDist);
    for j = 1:size(G,1)
        theta = 180*G(j,3)/pi;
        theta = theta + 90; % range: -90 to 90 is transferred to 0 to 180
        index = ceil(theta/binDist);
        if (index <= 0)
            index = 1;
        elseif (index > 180/binDist)
            index = 180/binDist;
        end
        H(1,index) = H(1,index)+1;
    end
    
function Hind = findIndividualHistogram(G,binDist)
i = 1;
Hind = [];
while (i < size(G,1))
    curve_num = G(i,4);
    idx = find(G(:,4) == curve_num);
    g = G(idx,:);
    H = findHistogram(g,binDist);
    Hind = [Hind;H];
    i = i+size(g,1);
end
    %}
function [hThresh Hhist] = find_histogram(Hue);
hsSize = 0.1;
multi = 10;
mn = multi*min(Hue);
mx = multi*max(Hue);
minH = floor(mn)/multi;
maxH = ceil(mx)/multi;
hThresh = [minH:hsSize:maxH];
nBin = size(hThresh,2)-1;
if (hThresh(1,nBin+1) < maxH)
    nBin = nBin+1;
    hThresh(1,nBin+1) = hThresh(1,nBin)+hsSize;
    maxH = hThresh(1,nBin+1);
end
Hhist = zeros(1,nBin);

for i = 1:nBin
    Hhist(1,i) = sum(Hue >= hThresh(1,i) & Hue < hThresh(1,i+1));
end

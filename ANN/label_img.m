function Output = label_img(gt, Data, numIdx)
% label Data using the barycenter of superpatch of gt

Y = regionprops(Data, 'centroid');
barycenter = round(cat(1, Y.Centroid));
Output = zeros(numIdx, 1);

for i = 1:numIdx
   Output(i) = gt(barycenter(i, 2), barycenter(i, 1)); 
end

end


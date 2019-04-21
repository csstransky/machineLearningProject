function SegGrad = calculate_segmentGradient(curve, extremum);
SegGrad = [];
for i = 1:size(extremum,2)-1
    SegGrad = [SegGrad atan((curve(extremum(i+1),2) - curve(extremum(i),2))/(curve(extremum(i+1),1) - curve(extremum(i),1)))];
end
 
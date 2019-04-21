function [orinetation edgeimage] = write_orientation_images(curve_sm, gard_curve, curve_num, sizex, sizey)
%fid = fopen('gradient_524000e_5251000n_LB_02.txt','w');

orinetation = zeros(sizex,sizey);
edgeimage = zeros(sizex,sizey);

for i = 1:curve_num
    Pts = round(curve_sm{i});
    x = sizex/2 - Pts(:,2);
    y = sizey/2 + Pts(:,1);
    for j = 1:size(Pts,1)
        orinetation(x(j,1),y(j,1)) = gard_curve{i}(j,1);
        edgeimage(x(j,1),y(j,1)) = i;
    end
end
%figure; imshow(mat2gray(orinetation));
%figure; imshow(mat2gray(edgeimage));
function [orinetation edgeimage] = extractImageLines(A,R,useAdjustment,useRefinement,fnameout)
% add paths to the Matlab path list
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\Digital Ortho Photo');
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\LAS Files');
%path(path,'C:\Develop\Knoch');
%path(path,'C:\Awrangjeb\Drive I\Hobert');

% read orhto image
%[A, R, bbox] = geotiffread('524000e_5251000n_LB_02.tif');

%extra = 0;
%ImgRect = [bbox(1,1)-extra,bbox(1,2)-extra; bbox(1,1)-extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(1,2)-extra];
%extra = 25;
%ImgRect1 = [bbox(1,1)-extra,bbox(1,2)-extra; bbox(1,1)-extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(1,2)-extra];

%mapshow(A, R); hold on; plot3(ALS(:,1), ALS(:,2), ALS(:,3),'.');
%Img = read_imagePixels('image_decor_crop7.txt');
%Img = Img';
%A = decorrstretch(A1,'Tol',0.01);% figure; imshow(C)
%mask = imread('mask_ground_e3490n58055_TL_sample_dem1_pixels.jpg');

Ig = rgb2gray(A(:,:,1:3));
Res = abs(R(2,1));
minLen = 3/Res;

%fp = fopen('mask_ground_crop7.txt', 'r');
%size = fscanf(fp, '%d%d%d',[1 3]);
%maskData = fscanf(fp, '%f %f %f', [3 inf]);
%fclose(fp);
%maskData = maskData';

%read ALS data
%fp = fopen('last_ground_crop7.txt', 'r');
%ALS = fscanf(fp, '%f %f %f', [3 inf]);
%fclose(fp);
%ALS = ALS';

%read Sigma(NDVI)
%fp = fopen('image_ndvi_sigma_crop7.txt', 'r');
%sz = fscanf(fp, '%d%d %d',[1 3]);
%N = fscanf(fp, '%f %f %f', [3 inf]);
%fclose(fp);
%N = N';

% convert image to gray
%I = rgb2gray(A);
%imshow(A);figure; imshow(I);
% Detect edges
%BW = edge(I,'canny');  
%[ght mask] = find_gHeightMask2(ALS,A,R,N);
%[ght mask ALSa maskg ALSg] = find_gHeightMask4(ALS,A,R);
%mask = maskg; ALSa = ALSg;
BW = edge(Ig,'canny'); 

% Extract curves from the edge-image
Gap_size = 1; minLength = minLen;
%[BW_sm curve_sm] = smooth_edgeImage01(BW,Gap_size,minLength);
[curve,curve_start,curve_end,curve_mode,curve_num,TJ,img1] = extract_curve(BW, Gap_size, minLength);  

%figure;
%imshow(img1); hold on;
%mapshow(img1,R); 
%hold on; plot3(ALS(:,1), ALS(:,2),ALS(:,3),'.');

%find smoothed curves and calculate their tangent(gradiants) on each point
[sizex sizey sizez] = size(A);
[smoothed_curve gard_curve c_intersect] = smooth_gradient01(curve, curve_mode, curve_num, sizex, sizey);
%[smoothed_curve grad_curve c_intersect] = smooth_gradient(curve, curve_mode, curve_num, R);

if useRefinement
    [orinetation edgeimage] = write_orientation_images(smoothed_curve, gard_curve, curve_num, sizex, sizey);
else
    orinetation = [];
    edgeimage = [];
end

if useAdjustment
%
minSegLength = minLen;
%gradThresh = 0.1;
%figure; imshow(~BW_sm); hold on;
%NgbrDiff = 0.1; Exclude = 5;
for i=1:curve_num    
    extremum{i} = selectPoint(curve{i}, curve_mode(i,1), curve_num, sizex, sizey);
    %showExtremum(curve{i},extremum{i});
    SegGrad{i} = calculate_segmentGradient(smoothed_curve{i}, extremum{i}); % not necessary
    x = smoothed_curve{i}(:,1);
    y = smoothed_curve{i}(:,2);
    %grad = gard_curve{i};
    L = size(x,1);
    %j = 1;
    %i1 = i
    Info{i} = []; %ind = 1;
    %creatNew = 1; tag = 0;    
    for j=1:size(SegGrad{i},2)
        %grad_diff(i,j) = abs(mean(grad(extremum{i}(1,j):extremum{i}(1,j+1)) - SegGrad{i}(1,j)));
        %plot([curve{i}(extremum{i}(1,j),2) curve{i}(extremum{i}(1,j+1),2)], [curve{i}(extremum{i}(1,j),1) curve{i}(extremum{i}(1,j+1),1)], ':r');
        if ((extremum{i}(1,j+1)-extremum{i}(1,j)) > minSegLength)
            Info{i}(j,1) = extremum{i}(1,j);
            Info{i}(j,2) = extremum{i}(1,j+1);
            %plot([curve{i}(extremum{i}(1,j),2) curve{i}(extremum{i}(1,j+1),2)], [curve{i}(extremum{i}(1,j),1) curve{i}(extremum{i}(1,j+1),1)], '-m','linewidth',2);
        else
            %Info{i}(j,1:2) = [];
        end
    end 
end
%hold off;

%show straight line segments
%figure; imshow(img1); hold on;
%for i=1:curve_num
%   mat = Info{i};
%    if (size(mat,1)>0)
%        for j=1:size(mat,1)
%            if (mat(j,1) > 0)
            %plot(curve{i}(mat(j,1),2),curve{i}(mat(j,1),1),'+r');
            %plot(curve{i}(mat(j,2),2),curve{i}(mat(j,2),1),'+r'); hold on;
%            plot([curve{i}(mat(j,1),2) curve{i}(mat(j,2),2)], [curve{i}(mat(j,1),1) curve{i}(mat(j,2),1)],'-r'); hold on;
            %out = [i j curve{i}(mat(j,1),:) curve{i}(mat(j,2),:)]
%            end
%       end
%    end
%end
%hold off;


% Straight-line fit
%msk = ones(size(maskg));
%figure; mapshow(msk,R); hold on;
%figure; mapshow(Ig,R); hold on;
fp = fopen(fnameout,'w');
for i=1:curve_num
    x = sizex/2 - smoothed_curve{i}(:,2);
    y = sizey/2 + smoothed_curve{i}(:,1);     
    mat = Info{i};
    O = ones(size(x,1),1);
    if (size(mat,1)>0)
        for j=1:size(mat,1)
            if (mat(j,1) > 0)
                xs = x(mat(j,1):mat(j,2),1);
                ys = y(mat(j,1):mat(j,2),1);
                Os = O(mat(j,1):mat(j,2),1);
                P = [xs ys Os]*R;
                Ls = size(P,1);
                Coef = polyfit(P(:,1),P(:,2),1);
                
                m = Coef(1,1);
                c = Coef(1,2);
                if (c == Inf || c == -Inf)
                    here = 1;
                else
                    mp = -1/m;
                    c1 = P(1,2)-mp*P(1,1);
                    c2 = P(Ls,2)-mp*P(Ls,1);

                    C1(1,1) = (c1-c)/(m-mp);
                    C1(1,2) = m*C1(1,1)+c;

                    C2(1,1) = (c2-c)/(m-mp);
                    C2(1,2) = m*C2(1,1)+c;

                     %plot([C1(1,1) C2(1,1)], [C1(1,2) C2(1,2)], '-g','linewidth',2); hold on;

                     y1 = floor((C1(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C1(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                     x1 = floor((C1(1,2) - R(3,2) - y1*R(2,2))/R(1,2));

                     y2 = floor((C2(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C2(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                     x2 = floor((C2(1,2) - R(3,2) - y2*R(2,2))/R(1,2));
                     %plot([y1 y2], [x1 x2], '-g','linewidth',2); hold on;

                     Info{i}(j,11:12) = C1;
                     Info{i}(j,13:14) = C2;
                     Info{i}(j,15:16) = [m c];
                     fprintf(fp,'%f\t %f\t %f\t %f\t %f\t %f\n',C1(1,1), C1(1,2), C2(1,1),C2(1,2),m,c);
                end
                
            end
        end
    end
end
fclose(fp);
%hold off;
%
end%if useAdjustment
here = 1;
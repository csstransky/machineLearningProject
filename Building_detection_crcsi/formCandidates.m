function formCandidates(A, R, bbox, maskg, Rm, ndvi, ndviThresh, eMask, useEntropy, fnamein, fnameout)
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\Digital Ortho Photo');
%path(path,'C:\Awrangjeb\Drive I\Sydney-ALS\LAS Files');
%path(path,'C:\Develop\Knoch');
%path(path,'C:\Awrangjeb\Drive I\Hobert');

% read orhto image
%[A, R, bbox] = geotiffread('524000e_5251000n_LB_02.tif');
extra = 0;
ImgRect = [bbox(1,1)-extra,bbox(1,2)-extra; bbox(1,1)-extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(1,2)-extra];
%extra = 25;
%ImgRect1 = [bbox(1,1)-extra,bbox(1,2)-extra; bbox(1,1)-extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(2,2)+extra; bbox(2,1)+extra,bbox(1,2)-extra];

%read Line data
fp = fopen(fnamein, 'r');
L = fscanf(fp, '%f %f %f %f %f %f', [12 inf]);
fclose(fp);
L = L';

%maskg = imread('mask_ground_524000e_5251000n_LB_02_dem1_pixels.jpg');
%figure; mapshow(maskg,R); hold on;
%for i=1:size(L,1)
%    P1 = L(i,1:2);
%    P2 = L(i,3:4);
%    plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)], '-c');
    
%    CP =  [(P1(1,1) + P2(1,1))/2 (P1(1,2) + P2(1,2))/2];
%    text(CP(1,1),CP(1,2),int2str(i), 'Color', 'cyan');
%end
%hold off;

%{
%read Mask data
fp = fopen('mask_ground_524000e_5251000n_LB_02_dem1_pixels.txt', 'r');
size1 = fscanf(fp, '%f %f %f', [1 3]);
M = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
M = M';

%read NDVI_sigma data
fp = fopen('image_ndvi_524000e_5251000n_LB_02.txt', 'r');
size2 = fscanf(fp, '%f %f %f', [1 3]);
N = fscanf(fp, '%f %f %f', [3 inf]);
fclose(fp);
N = N';

%read texture info
fp = fopen('image_texture_524000e_5251000n_LB_02.txt','r');
T = fscanf(fp, '%f %f %d %d', [4 inf]);
fclose(fp);
T = T';

%read mask
mask = imread('mask_up_524000e_5251000n_LB_02_dem1_pixels.jpg');
%}

% filter using ALS
dAl = 1.5; % radius of circle
%htThresh = 9.5;%7.5;
edgeThresh = 0.70; %nBin = 16; 
%Sin = 0.50; % Spacing_in_flight
%Sacross = 0.50; % Spacing_across_flight
%figure; mapshow(A,R); hold on; 
%ndviThresh = 10; 
entropyThresh = 0.30; 
%edgeShadowThresh = 0.05;

[srtLen idx] = sort(L(:,7),'descend');
 %seg_num = 0; Segs = []; dAc = 6; dAl = 3;

seg_num = size(L,1);
%max_buildingLength = 50;
min_buildingLength = 3;
%half_diagonSqr = max_buildingLength*max_buildingLength/2;
Buildings = []; bNumber = 0; %mThresh = 0.05;
Tag = zeros(seg_num,1);

%figure; mapshow(A(:,:,1:3),R); hold on;
%figure; mapshow(A(:,:,1:3),R); hold on;
for t = 1:seg_num    
    i = idx(t,1); % choose long lines first
    %if (i == 269)
    %    here = 1;
    %end
    if (Tag(i,1) < 1) %if it is not used yet
        %plot([Segs(i,1) Segs(i,3)], [Segs(i,2) Segs(i,4)], '-g'); hold on;
        %plot(Segs(i,11), Segs(i,12),'+g');
         %plot([L(i,1) L(i,3)], [L(i,2) L(i,4)], '-g', 'linewidth',2); hold on;
         %plot(L(i,8),L(i,9),'*b');
        ret = testInBuidlings_01(Buildings,L(i,:));
        %i1 = i
        %if (i == 95)
        %    here = 1;
        %end
        if (ret)
            bNumber = bNumber+1;        
            Tag(i,1) = bNumber;
            Seg1 = L(i,:);
            FlagT = L(i,12);
            P1 = Seg1(1,1:2);
                P2 = Seg1(1,3:4);
                
                %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-k','LineWidth', 2); hold on;
                %plot(Seg1(1,8),Seg1(1,9),'og'); hold on;
                
                m12 = Seg1(1,5);
                c12 = Seg1(1,6);
                % process L34
                [Q1 Q2 Q3 Q4] = find_rectVertices(m12,c12,dAl,P1,P2);
                Y = m12*[Seg1(1,8); Q1(1,1)] + c12;
                Y_diff = Y-[Seg1(1,9); Q1(1,2)];
                if ((Y_diff(1,1) <= 0 && Y_diff(2,1) <= 0) || (Y_diff(1,1)  > 0 && Y_diff(2,1) > 0)) % both C1 of Seg2 and in-point of Seg1 are in same side
                    Q1in = Q1;
                    Q3in = Q3;
                    
                    Q2in = Q2;
                    Q4in = Q4;                    
                else
                    Q1in = Q2;
                    Q3in = Q4;
                    
                    Q2in = Q1;
                    Q4in = Q3;                    
                end
                
                [Q1in Q3in] = extend_initial_building(A,R,maskg,Rm,-1/m12,dAl,Q1in,Q3in,P1,ndvi,eMask,useEntropy,edgeThresh,ndviThresh, entropyThresh, ImgRect,FlagT); % extend left
                P3 = Q3in; P4 = Q1in;                               
                [P2 P3] = extend_initial_building(A,R,maskg,Rm,m12,dAl,P2,P3,P1,ndvi,eMask,useEntropy,edgeThresh,ndviThresh, entropyThresh, ImgRect,FlagT); % extend left                
                [P4 P1] = extend_initial_building(A,R,maskg,Rm,m12,dAl,P4,P1,P2,ndvi,eMask,useEntropy,edgeThresh,ndviThresh, entropyThresh, ImgRect,FlagT); % extend left
                
                 
                 Len1 = ceil(sqrt((P1(1,1)-P4(1,1))*(P1(1,1)-P4(1,1)) + (P1(1,2)-P4(1,2))*(P1(1,2)-P4(1,2))));
                Len2 = ceil(sqrt((P1(1,1)-P2(1,1))*(P1(1,1)-P2(1,1)) + (P1(1,2)-P2(1,2))*(P1(1,2)-P2(1,2))));
                if (Len1 < min_buildingLength || Len2 < min_buildingLength)
                    %discard Seg1 as building edge
                    Tag(i,1) = 9999;                    
                    bNumber = bNumber-1;
                    
                    %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-r','LineWidth', 3); hold on;
                    %plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-r','LineWidth', 3); hold on;
                    %plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-r','LineWidth', 3); hold on;
                    %plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-r','LineWidth', 3); hold on;   
                else
                    %ret1 = testInRectanglesAxisParallel(ImgRect1(1,:), ImgRect1(2,:), ImgRect1(3,:), ImgRect1(4,:), [P1;P2;P3;P4]);
                    %if (ret1)
                        %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-b','LineWidth', 3); hold on;
                        %plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-b','LineWidth', 3); hold on;
                        %plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-b','LineWidth', 3); hold on;
                        %plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-b','LineWidth', 3); hold on;   
                    %end
                    Buildings = [Buildings;[P1 P2 P3 P4 FlagT]];
                end                      
        else %if(ret)
            Tag(i,1) = 9999; % for unused line segment
        end%if(ret)
    end%if
    %out = t
end%for
%hold off;

fp = fopen(fnameout,'w');
for i = 1:size(Buildings,1)
    fprintf(fp,'%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %d\n', Buildings(i,1), Buildings(i,2), Buildings(i,3), Buildings(i,4), Buildings(i,5), Buildings(i,6), Buildings(i,7), Buildings(i,8), Buildings(i,9));
end
fclose(fp);
here = 1;
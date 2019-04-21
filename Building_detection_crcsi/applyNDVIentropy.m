function applyNDVIentropy(A, R, mask, Rm, ndvi, ndviThresh, eMask, useEntropy, useAdjustment, fnamein, fnameout)

%{
%compute NDVI
[ndvi ndvisig] = computeNDVI(A, 0);

[height width bands] = size(A);
if bands == 3 %we need sigma of pseudo NDVI, so convert ndvisig to 8-bit greyscale image
    mn = min(min(ndvisig));
    mx = max(max(ndvisig));
    ndvisig = uint8((ndvisig-mn)*255/(mx-mn)+0.5);    
    ndviThresh = 48;
    nvdi = ndvisig;
else
    ndviThresh = 10;
end


if useEntropy  
    % calculate image texture
    Ig = rgb2gray(A(:,:,1:3));
    E = entropyfilt(Ig);
    Eim = mat2gray(E);    
    eMask = im2bw(Eim, .8);
    %figure; imshow(eMask);
end
%}

%read Line data
fp2 = fopen(fnamein, 'r');
L = fscanf(fp2, '%f %f %f %f %f %f', [6 inf]);
fclose(fp2);
L = L';

%first-last difference
%fp = fopen('first_last_diff_crop7.txt', 'r');
%Ad = fscanf(fp, '%f %f %f', [3 inf]);
%fclose(fp);
%Ad = Ad';
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
%}

%read mask
%mask = imread('mask_ground_524000e_5251000n_LB_02_dem1_pixels.jpg');

% filter using ALS
dAl = 2; % radius of circle 
%dAl1 = 5; dThresh = 2;
%htThresh = 9.5;%7.5;
edgeThresh = 0.70; %nBin = 16; 
%Sin = 0.50; % Spacing_in_flight
%Sacross = 0.50; % Spacing_across_flight

%ndvi = imread('e3490n58055_TL_sigNDVI8bit.tif');
%figure; mapshow(A(:,:,1:3),R); hold on; 
%figure; imshow(A(:,:,1:3)); hold on; 
L1 = []; L2 = []; 

entropyThresh = 0.30; %edgeShadowThresh = 0.05;
%fp1 = fopen('TextureLines.txt','w'); 
num = 0;
for i = 1:size(L,1)
%for i = 1:50
    %mat = Info{i};
    %if (size(mat,1)>0)
        %for j=1:size(mat,1)
            %if (mat(j,1) > 0) % if at least one segment is found               
                 C1 = L(i,1:2);
                 C2 = L(i,3:4);
                 m = L(i,5);
                 c = L(i,6);
                  
                 [P1 P2 P3 P4] = find_rectVertices(m,c,dAl,C1,C2);
                 %[P1a P2a P3a P4a] = find_rectVertices(m,c,dAl1,C1,C2);
                
                %plot([C1(1,1) C2(1,1)], [C1(1,2) C2(1,2)],'-g'); hold on;
                %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-k'); hold on;
                %plot([P2(1,1) P4(1,1)], [P2(1,2) P4(1,2)],'-r'); hold on;
                %plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-k'); hold on;
                %plot([P3(1,1) P1(1,1)], [P3(1,2) P1(1,2)],'-b'); hold on;
                %plot([P1a(1,1) P2a(1,1)], [P1a(1,2) P2a(1,2)],'-c'); hold on;
                %plot([P2a(1,1) P4a(1,1)], [P2a(1,2) P4a(1,2)],'-c'); hold on;
                %plot([P3a(1,1) P4a(1,1)], [P3a(1,2) P4a(1,2)],'-c'); hold on;
                %plot([P3a(1,1) P1a(1,1)], [P3a(1,2) P1a(1,2)],'-c'); hold on;
                
                 %out = [i j]
                 
                 %convert geographic to image
                    %for C1
                    Pout = obj2obs(C1,Rm);
                    c1x = Pout(1,1); c1y = Pout(1,2);
                 %c1y = floor((C1(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C1(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                 %c1x = floor((C1(1,2) - R(3,2) - c1y*R(2,2))/R(1,2));
                    %for C2
                    Pout = obj2obs(C2,Rm);
                    c2x = Pout(1,1); c2y = Pout(1,2);
                 %c2y = floor((C2(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C2(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                 %c2x = floor((C2(1,2) - R(3,2) - c2y*R(2,2))/R(1,2));
                    %for P1
                 %p1y = floor((P1(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - P1(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                 %p1x = floor((P1(1,2) - R(3,2) - c1y*R(2,2))/R(1,2));
                 
                 len = sqrt((c1x-c2x)*(c1x-c2x) + (c1y-c2y)*(c1y-c2y));
                 %wid = sqrt((c1x-p1x)*(c1x-p1x) + (c1y-p1y)*(c1y-p1y));
                 %nPoints or area
                 %nPoints = len*wid;
                 %Info{i}(j,17) = len;
                 %Info{i}(j,18) = nPoints;
                 
                %[P1 P2 P3 P4] = find_rectVertices(m,c,dAl,C1,C2);
                
                %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-g'); hold on;
                %plot([P2(1,1) P4(1,1)], [P2(1,2) P4(1,2)],'-b'); hold on;
                %plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-y'); hold on;
                %plot([P3(1,1) P1(1,1)], [P3(1,2) P1(1,2)],'-r'); hold on;
               
                %plot([P1(1,1) C1(1,1)], [P1(1,2) C1(1,2)],':b','linewidth',2); hold on;
                %plot([P1(1,1) P3(1,1)], [P1(1,2) P3(1,2)],':b','linewidth',2); hold on;
                %plot([P3(1,1) C2(1,1)], [P3(1,2) C2(1,2)],':b','linewidth',2); hold on;
                
                %plot([P2(1,1) C1(1,1)], [P2(1,2) C1(1,2)],':c','linewidth',2); hold on;
                %plot([P2(1,1) P4(1,1)], [P2(1,2) P4(1,2)],':c','linewidth',2); hold on;
                %plot([P4(1,1) C2(1,1)], [P4(1,2) C2(1,2)],':c','linewidth',2); hold on;
                
                C1m = obj2obs(C1,Rm);
                C2m = obj2obs(C2,Rm);
                P1m = obj2obs(P1,Rm);
                P2m = obj2obs(P2,Rm);
                P3m = obj2obs(P3,Rm);
                P4m = obj2obs(P4,Rm);
                
                %plot([C1m(1,2) C2m(1,2)], [C1m(1,1) C2m(1,1)],'-g'); hold on;
                %plot([P1m(1,2) P2m(1,2)], [P1m(1,1) P2m(1,1)],'-k'); hold on;
                %plot([P2m(1,2) P4m(1,2)], [P2m(1,1) P4m(1,1)],'-r'); hold on;
                %plot([P3m(1,2) P4m(1,2)], [P3m(1,1) P4m(1,1)],'-k'); hold on;
                %plot([P3m(1,2) P1m(1,2)], [P3m(1,1) P1m(1,1)],'-b'); hold on;
                
                [totalaP numValaP] = findImagePointCounts(mask, P1m, C1m, C2m, P3m, 0);
                [totalbP numValbP] = findImagePointCounts(mask, C1m, P2m, P4m, C2m, 0);
                
                %M_above = find_ALS(A,R,M,P1,C1,C2,P3);
                %M_below = find_ALS(A,R,M,C1,P2,P4,C2);
                
                %above = sum((M_above(:,3) == 0))/nPoints;
                %below = sum((M_below(:,3) == 0))/nPoints;              
                nPoints = max([totalaP, totalbP]);
                above = numValaP/nPoints;
                below = numValbP/nPoints;              
                
                %mpx = (C1(1,1) + C2(1,1))/2;
                %mpy = (C1(1,2) + C2(1,2))/2;
                
                %plot([C1(1,1) C2(1,1)], [C1(1,2) C2(1,2)],'-k','linewidth',2); hold on;
                %text(mpx, mpy, int2str(i),'Color','b');
                %if (i == 11)
                %   here = 1;
                %end
                
                if (above < edgeThresh && below < edgeThresh)
                    % line segment is NOT an edge of roof and avoid it
                    
                elseif (above >= edgeThresh)
                    
                    C1i = obj2obs(C1,R);
                    C2i = obj2obs(C2,R);
                    P1i = obj2obs(P1,R);
                    %P2i = obj2obs(P2,R);
                    P3i = obj2obs(P3,R);
                    %P4i = obj2obs(P4,R);
                
                    %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-g'); hold on;
                    %plot([P1i(1,2) P2i(1,2)], [P1i(1,1) P2i(1,1)],'-k'); hold on;
                    %plot([P2i(1,2) P4i(1,2)], [P2i(1,1) P4i(1,1)],'-r'); hold on;
                    %plot([P3i(1,2) P4i(1,2)], [P3i(1,1) P4i(1,1)],'-k'); hold on;
                    %plot([P3i(1,2) P1i(1,2)], [P3i(1,1) P1i(1,1)],'-b'); hold on;
                
                    m = findImagePointMean(ndvi, P1i, C1i, C2i, P3i);
                    
                    % the line segment is the roof-edge
                    %N_above = find_ALS(A,R,N,P1,C1,C2,P3);
                    %D_above = find_ALS(A,R,Ad,P1,P1a,P3a,P3);
                    %m = mean(N_above(:,3));
                    %Dinfo = sum(D_above(:,3)>= dThresh)/size(D_above,1);
                    if (m <= ndviThresh)
                        mpxx = (P1(1,1) + P3(1,1))/2;
                        mpyy = (P1(1,2) + P3(1,2))/2;
                        L1 = [L1; [L(i,1:6) len mpxx mpyy 0]];
                        %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-c','linewidth',2); hold on;
                        %text(mpx, mpy, int2str(L(i,7)),'Color','c');
                        %plot(mpx,mpy,'ok', 'LineWidth',1, 'MarkerEdgeColor','b', 'MarkerFaceColor','b', 'MarkerSize',6); hold on;
                    elseif useEntropy
                        m1 = (P3(1,2)-P1(1,2))/(P3(1,1)-P1(1,1));
                        c1 = P1(1,2)-m1*P1(1,1);
                        [P1a P2a P3a P4a] = find_rectVertices(m1,c1,1,P1,P3);
                        
                        P1ai = obj2obs(P1a,R);
                        P2ai = obj2obs(P2a,R);
                        P3ai = obj2obs(P3a,R);
                        P4ai = obj2obs(P4a,R);
                        
                        %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-g'); hold on;
                        %plot([P1ai(1,2) P2ai(1,2)], [P1ai(1,1) P2ai(1,1)],'-k'); hold on;
                        %plot([P2ai(1,2) P4ai(1,2)], [P2ai(1,1) P4ai(1,1)],'-r'); hold on;
                        %plot([P3ai(1,2) P4ai(1,2)], [P3ai(1,1) P4ai(1,1)],'-k'); hold on;
                        %plot([P3ai(1,2) P1ai(1,2)], [P3ai(1,1) P1ai(1,1)],'-b'); hold on;
                        
                        [totaleP numValeP] = findImagePointCounts(eMask, P1ai, P3ai, P4ai, P2ai, 1);
                        m_entr = numValeP/totaleP;
                        
                        %{
                        %plot([P1a(1,1) P2a(1,1)], [P1a(1,2) P2a(1,2)],'-c'); hold on;
                        %plot([P2a(1,1) P4a(1,1)], [P2a(1,2) P4a(1,2)],'-c'); hold on;
                        %plot([P3a(1,1) P4a(1,1)], [P3a(1,2) P4a(1,2)],'-c'); hold on;
                        %plot([P3a(1,1) P1a(1,1)], [P3a(1,2) P1a(1,2)],'-c'); hold on;
                        M_abv = find_ALS(A,R,M,P1a,P3a,P4a,P2a);
                        abv = sum((M_abv(:,3) == 0))/size(M_abv,1);
                        N_abv = find_ALS(A,R,N,P1a,P3a,P4a,P2a);
                        ma = mean(N_abv(:,3));
                        T_abv = find_ALS(A,R,T,P1a,P3a,P4a,P2a);
                        id_edge = (T_abv(:,3) == 1);
                        id_entr = (T_abv(:,4) == 1);
                        m_edge = sum(id_edge)/size(T_abv,1);
                        m_entr = sum(id_entr)/size(T_abv,1);
                        
                        T_above = find_ALS(A,R,T,P1,C1,C2,P3);
                        id_edges = (T_above(:,3) == 1);
                        id_entro = (T_above(:,4) == 1);
                        m_edges = sum(id_edges)/size(T_above,1);
                        m_entro = sum(id_entro)/size(T_above,1);
                        %}
                        num = num+1;
                        
                        %if (abv >= above && m_entr <= entropyThresh && m_edge > edgeShadowThresh)
                        if (m_entr <= entropyThresh)
                            mpxx = (P1(1,1) + P3(1,1))/2;
                            mpyy = (P1(1,2) + P3(1,2))/2;
                            L2 = [L2; [L(i,1:6) len mpxx mpyy 1]];                        
                            %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-c','linewidth',2); hold on;
                             %text(mpx, mpy, int2str(L(i,7)),'Color','c');
                            %text(mpx,mpy,int2str(num),'Color','k'); hold on;                                             
                            %fprintf(fp1,'%d\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n', num, above, abv, m, ma, m_edges, m_edge, m_entro, m_entr);
                        end
                        %plot(mpx,mpy,'ok', 'LineWidth',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'MarkerSize',6); hold on;
                        %plot(mpx,mpy,'ok', 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4); hold on;
                    end                    
                else
                    
                    C1i = obj2obs(C1,R);
                    C2i = obj2obs(C2,R);
                    %P1i = obj2obs(P1,R);
                    P2i = obj2obs(P2,R);
                    %P3i = obj2obs(P3,R);
                    P4i = obj2obs(P4,R);
                
                    %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-g'); hold on;
                    %plot([P1i(1,2) P2i(1,2)], [P1i(1,1) P2i(1,1)],'-k'); hold on;
                    %plot([P2i(1,2) P4i(1,2)], [P2i(1,1) P4i(1,1)],'-r'); hold on;
                    %plot([P3i(1,2) P4i(1,2)], [P3i(1,1) P4i(1,1)],'-k'); hold on;
                    %plot([P3i(1,2) P1i(1,2)], [P3i(1,1) P1i(1,1)],'-b'); hold on;
                    
                    m = findImagePointMean(ndvi, C1i, P2i, P4i, C2i);
                    
                    %N_below = find_ALS(A,R,N,C1,P2,P4,C2);
                    %D_below = find_ALS(A,R,Ad,P2a,P2,P4,P4a);
                    %m = mean(N_below(:,3));                    
                    %Dinfo = sum(D_below(:,3)>= dThresh)/size(D_below,1);
                    if (m <= ndviThresh)
                        mpxx = (P2(1,1) + P4(1,1))/2;
                        mpyy = (P2(1,2) + P4(1,2))/2;
                        L1 = [L1; [L(i,1:6) len mpxx mpyy 0]];
                        %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-c','linewidth',2); hold on;                       
                        %text(mpx, mpy, int2str(L(i,7)),'Color','c');
                        %plot(mpx,mpy,'ok', 'LineWidth',1, 'MarkerEdgeColor','b', 'MarkerFaceColor','b', 'MarkerSize',6); hold on;
                    elseif useEntropy
                        m1 = (P4(1,2)-P2(1,2))/(P4(1,1)-P2(1,1));
                        c1 = P2(1,2)-m1*P2(1,1);
                        [P1a P2a P3a P4a] = find_rectVertices(m1,c1,1,P2,P4);
                        
                        P1ai = obj2obs(P1a,R);
                        P2ai = obj2obs(P2a,R);
                        P3ai = obj2obs(P3a,R);
                        P4ai = obj2obs(P4a,R);
                        
                        %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-g'); hold on;
                        %plot([P1ai(1,2) P2ai(1,2)], [P1ai(1,1) P2ai(1,1)],'-k'); hold on;
                        %plot([P2ai(1,2) P4ai(1,2)], [P2ai(1,1) P4ai(1,1)],'-r'); hold on;
                        %plot([P3ai(1,2) P4ai(1,2)], [P3ai(1,1) P4ai(1,1)],'-k'); hold on;
                        %plot([P3ai(1,2) P1ai(1,2)], [P3ai(1,1) P1ai(1,1)],'-b'); hold on;
                        
                        [totaleP numValeP] = findImagePointCounts(eMask, P1ai, P3ai, P4ai, P2ai, 1);
                        m_entr = numValeP/totaleP;
                        
                        %{
                        %plot([P1a(1,1) P2a(1,1)], [P1a(1,2) P2a(1,2)],'-c'); hold on;
                        %plot([P2a(1,1) P4a(1,1)], [P2a(1,2) P4a(1,2)],'-c'); hold on;
                        %plot([P3a(1,1) P4a(1,1)], [P3a(1,2) P4a(1,2)],'-c'); hold on;
                        %plot([P3a(1,1) P1a(1,1)], [P3a(1,2) P1a(1,2)],'-c'); hold on;
                        M_blw = find_ALS(A,R,M,P1a,P3a,P4a,P2a);
                        blw = sum((M_blw(:,3) == 0))/size(M_blw,1);
                        N_blw = find_ALS(A,R,N,P1a,P3a,P4a,P2a);
                        ma = mean(N_blw(:,3));
                        T_blw = find_ALS(A,R,T,P1a,P3a,P4a,P2a);
                        id_edge = (T_blw(:,3) == 1);
                        id_entr = (T_blw(:,4) == 1);
                        m_edge = sum(id_edge)/size(T_blw,1);
                        m_entr = sum(id_entr)/size(T_blw,1);
                        
                        T_below = find_ALS(A,R,T,C1,P2,P4,C2);
                        id_edges = (T_below(:,3) == 1);
                        id_entro = (T_below(:,4) == 1);
                        m_edges = sum(id_edges)/size(T_below,1);
                        m_entro = sum(id_entro)/size(T_below,1);                        
                        %}
                        
                        num = num+1;
                        
                        %if (blw >= below && m_entr <= entropyThresh && m_edge > edgeShadowThresh)
                        if (m_entr <= entropyThresh)
                            mpxx = (P2(1,1) + P4(1,1))/2;
                            mpyy = (P2(1,2) + P4(1,2))/2;
                            L2 = [L2; [L(i,1:6) len mpxx mpyy 1]];
                            %plot([C1i(1,2) C2i(1,2)], [C1i(1,1) C2i(1,1)],'-c','linewidth',2); hold on;
                            %text(mpx, mpy, int2str(L(i,7)),'Color','c');
                            %text(mpx,mpy,int2str(num),'Color','k'); hold on;
                            %fprintf(fp1,'%d\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n', num, below, blw, m, ma, m_edges, m_edge, m_entro, m_entr);
                        end
                        %plot(mpx,mpy,'ok', 'LineWidth',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'MarkerSize',6); hold on;
                        %plot(mpx,mpy,'ok', 'LineWidth',2, 'MarkerEdgeColor','w', 'MarkerFaceColor','w', 'MarkerSize',8); hold on;
                        %plot(mpx,mpy,'ok', 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',4); hold on;
                    end
                end
            %else
             %   Info{i}(j,3) = 0; % line on roof
             %   Info{i}(j,4) = 0; % no significance                
            %end
        %end
    %end
    %out = [i m]
end
%hold off;
%fclose(fp1);

fp1 = fopen(fnameout,'w');
if useAdjustment
    %use image lines for adjustment
    %read image Line data
    fp = fopen('imageLines.txt', 'r');
    Li = fscanf(fp, '%f %f %f %f %f %f', [6 inf]);
    fclose(fp);
    Li = Li';
    
    [L1 F] = replaceImageLines(A,R,[L1;L2],Li);
    
    [L1 D V] = adjust_lines_image(A,R,L1,Li,F);    
    
    
    
    %compute priority
    maxLength = max(L1(:,7));
    %mxD = max(D);
    %L1(:,7) = L1(:,7)/maxLength + (mxD-D)/mxD + 3*F;
    L1(:,7) = L1(:,7)/maxLength + 3*F;
    
    
    %alternatively use local longest line for adjustment
    % uncomment the following line and comment all above code for
    % adjustment using image lines
    %[L1 F V] = adjust_lines(A,R,[L1;L2]); 
        
    for i = 1:size(L1,1)
        if (V(i,1) == 0)
            fprintf(fp1,'%f\t %f\t %f\t %f\t %f\t %f\t %f \t %f \t %f \t %d\n',L1(i,1), L1(i,2), L1(i,3),L1(i,4),L1(i,5), L1(i,6), L1(i,7),L1(i,8),L1(i,9), L1(i,10));
        end
    end
else
    L1 = [L1;L2];
    
    for i = 1:size(L1,1)
       fprintf(fp1,'%f\t %f\t %f\t %f\t %f\t %f\t %f \t %f \t %f \t %d\n',L1(i,1), L1(i,2), L1(i,3),L1(i,4),L1(i,5), L1(i,6), L1(i,7),L1(i,8),L1(i,9), L1(i,10));       
    end
end
fclose(fp1);

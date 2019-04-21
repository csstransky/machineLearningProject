function extractLines(mask, Rm, fname)

BW = edge(mask,'canny'); 

% Extract curves from the edge-image
Gap_size = 1; 
resolutionMask = 0.25;
minLength = 3/resolutionMask;

[BW_sm curve_sm] = smooth_edgeImage01(BW,Gap_size,minLength);
[curve,curve_start,curve_end,curve_mode,curve_num,TJ,img1] = extract_curve(BW_sm, Gap_size, minLength);  

%find smoothed curves and calculate their tangent(gradiants) on each point
[sizex sizey] = size(mask);
[smoothed_curve gard_curve c_intersect] = smooth_gradient01(curve, curve_mode, curve_num, sizex, sizey);

%
minSegLength = minLength;
for i=1:curve_num    
    extremum{i} = selectPoint(curve{i}, curve_mode(i,1), curve_num, sizex, sizey);
    SegGrad{i} = calculate_segmentGradient(smoothed_curve{i}, extremum{i}); % not necessary
    x = smoothed_curve{i}(:,1);
    y = smoothed_curve{i}(:,2);
    
    L = size(x,1);
    Info{i} = []; %ind = 1;
    for j=1:size(SegGrad{i},2)
        if ((extremum{i}(1,j+1)-extremum{i}(1,j)) > minSegLength)
            Info{i}(j,1) = extremum{i}(1,j);
            Info{i}(j,2) = extremum{i}(1,j+1);
            %plot([curve{i}(extremum{i}(1,j),2) curve{i}(extremum{i}(1,j+1),2)], [curve{i}(extremum{i}(1,j),1) curve{i}(extremum{i}(1,j+1),1)], '-m','linewidth',2);
        else
        end
    end 
end

%figure; mapshow(mask,Rm); hold on;
fp = fopen(fname,'w');
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
                P = [xs ys Os]*Rm;
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

                     %y1 = floor((C1(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C1(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                     %x1 = floor((C1(1,2) - R(3,2) - y1*R(2,2))/R(1,2));

                     %y2 = floor((C2(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C2(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                     %x2 = floor((C2(1,2) - R(3,2) - y2*R(2,2))/R(1,2));
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
here = 1;
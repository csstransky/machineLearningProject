function [L F] = replaceImageLines(A,R,L,Li)
%function adjust_lines();
%ResImage = abs(R(2,1));

%path(path,'I:\Sydney-ALS\Digital Ortho Photo');
%path(path,'I:/Sydney-ALS/ALS');
%path(path,'C:\Develop\Knoch');

% read orhto image
%[A, R, bbox] = geotiffread(Iname);

%read Line data
%fp = fopen('mask_ground_LSM_lines_mask_ndvi_crop6.txt', 'r');
%L = fscanf(fp, '%f %f %f %f %f %f', [9 inf]);
%fclose(fp);
%L = L';

%[sdist ids] = sort(L(:,7),'descend');
%L = L(ids,:);
%maxBuildingLength = 50;
%maxBuildingLengthSquare = maxBuildingLength *maxBuildingLength ;

for i = 1:size(L,1)
    C1 = L(i,1:2); C2 = L(i,3:4);
    Mp(i,:) = [(C1(1,1) + C2(1,1))/2 (C1(1,2) + C2(1,2))/2];
end

for i = 1:size(Li,1)
    C1 = Li(i,1:2); C2 = Li(i,3:4);
    Mpi(i,:) = [(C1(1,1) + C2(1,1))/2 (C1(1,2) + C2(1,2))/2];
end

%O = ones(size(Mpi,1),1);
%max_buildingLength = 50;
%max_sqr_dist = max_buildingLength * max_buildingLength;
F = zeros(size(L,1),1); % 0: never, 1: previously adjusted
%V = zeros(size(L,1),1); % 0: never, 1: previously adjusted
%D = -1*ones(size(L,1),1); % -1: not adjusted, +ve: distant to nearest parallel or perpendicular image line which was used to adjust
thetaThresh = 22.5;
%minLength = 60; %minimum length is 6m; 60 pixels for Knox
%read mask
%mask = imread('mask_ground_crop6.jpg');
disThersh = 3.0;
%figure; mapshow(A(:,:,1:3),R); hold on;
%figure; mapshow(mask,R); hold on;
%total = 0; 
T = size(L,1);
Ti = size(Li,1);

mi = Li(:,5);
ci = Li(:,6);
deni = sqrt((1+mi.*mi));
for i=1:T
    
    %if (total + sum(V) == T)
    %    break;
    %end
    %id = i;%find(L(:,11) == L(i,11)); % find lines with the same curve id       
    %L1 = L(i,:); % neigbor lines 
    %F1 = F(id,:); % neigbor flags
    %Mp1 = Mp(i,:); %neighbor midpoints
    %adsum = 0;
    %if (V(i,1) == 0 && L(i,7) >= minLength && F(i,1) == 0) %F(i,1) = 0: not visited yet; L(i,7) >= minLength: this line can be a ref line to fix others with the same curveID; V(i,1) = 0: this line is by default selected
    %if (V(i,1) == 0) %F(i,1) = 0: not visited yet; L(i,7) >= minLength: this line can be a ref line to fix others with the same curveID; V(i,1) = 0: this line is by default selected
        %difX = O*Mp1(1,1) - Mpi(:,1);
        %difY = O*Mp1(1,2) - Mpi(:,2);
        %cp_sqr = difX.*difX + difY.*difY;
        %[cpsqr idcps] = sort(cp_sqr,'ascend');
        
        d = abs(mi*Mp(i,1) - Mp(i,2) + ci)./deni;
        
        [ds ids] = sort(d,'ascend');
        
        %id = find(cp_sqr <= max_sqr_dist);
        %id = find(L(:,11) == L(i,11)); % find lines with the same curve id
  
        %L1 = L(id,:); % neigbor lines
        %F1 = F(id,:); % neigbor flags
        %Mp1 = Mp(id,:); %neighbor midpoints
        %num = size(F1,1) - sum(F1);
        %count = 0;        
        %[sdist sid] = sort(L1(:,7),'descend');
        %adsum = 0;
        %if (num > 0)
        %j = 1;
        %while(count < num-1 && j < size(L1,1)) % in F1: 0 = never, 1: previously adjusted, 2: newly adjusted
            %ind = sid(j,1); 
            
            %if (F1(j,1) ~= 2) %if not recently adjusted
                %m = L1(ind,5);
                %theta = 180*atan(m)/pi;
                %mper = -1/m;
                %thetaper = 180*atan(mper)/pi;
                
                %now = 1; adjust = 0;
                %while (adjust == 0 && cpsqr(now,1) <= maxBuildingLengthSquare && now <= size(Li,1) && now <= 10)%check first nearest lines only 
                %    idimg = idcps(now,1);
                
                
                P1 = L(i,1:2);
                P2 = L(i,3:4);                
                %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-b','linewidth',2); hold on;
                
                CP = Mp(i,1:2);
                %Pin = L1(1,8:9);
                P1c = [P1(1,1) - CP(1,1) P1(1,2) - CP(1,2)];
                P2c = [P2(1,1) - CP(1,1) P2(1,2) - CP(1,2)];
                %Pinc = [Pin(1,1) - CP(1,1) Pin(1,2) - CP(1,2)];
                
                m1 = (P2c(1,1) - P1c(1,1))/(P2c(1,2) - P1c(1,2));
                theta1 = 180*atan(m1)/pi;
                        
                k = 1; notFound = 1;
                while (k <= Ti && ds(k,1) < disThersh && notFound)
                %for k = 1:Ti               
                    %if (i ~= id(k,1) && F1(k,1) == 0)                     
                        idi = ids(k,1);
                        
                        Q1 = Li(idi,1:2);
                        Q2 = Li(idi,3:4); 
                        %plot([Q1(1,1) Q2(1,1)], [Q1(1,2) Q2(1,2)],'-m','linewidth',1); hold on;
                        
                        Q1c = [Q1(1,1) - CP(1,1) Q1(1,2) - CP(1,2)];
                        Q2c = [Q2(1,1) - CP(1,1) Q2(1,2) - CP(1,2)];
                        
                        m = (Q2c(1,1) - Q1c(1,1))/(Q2c(1,2) - Q1c(1,2));
                        theta = 180*atan(m)/pi;
                        %mper = -1/m;
                        %thetaper = 180*atan(mper)/pi;
                
                        %P1 = L1(k,1:2);
                        %P2 = L1(k,3:4); 
                        %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-m','linewidth',1); hold on;
                        %Pin = L1(k,8:9);
                        %P1c = [P1(1,1) - CP(1,1) P1(1,2) - CP(1,2)];
                        %P2c = [P2(1,1) - CP(1,1) P2(1,2) - CP(1,2)];
                        %Pinc = [Pin(1,1) - CP(1,1) Pin(1,2) - CP(1,2)];
                        
                        %m1 = (P2c(1,1) - P1c(1,1))/(P2c(1,2) - P1c(1,2));
                        %theta1 = 180*atan(m1)/pi;
                        thetadiff = abs(theta-theta1);
                        
                        
                        %Fl = 0;
                        if (thetadiff <= thetaThresh || abs(thetadiff-180) <= thetaThresh) % parallel to line at j of L1
                            %rot = theta-theta1;
                            %L1(k,:) = rotateLine(L1(k,:), rot);
                            %Fl = 1;
                            %adsum = adsum+1;
                        %else
                        %     thetadiff = abs(thetaper-theta1);
                        %     if (thetadiff <= thetaThresh || abs(thetadiff-180) <= thetaThresh) % perpendicular to line at j of L1
                        %            rot = thetaper-theta1;
                        %            %L1(k,:) = rotateLine(L1(k,:), rot);
                        %            Fl = 1;
                        %            adsum = adsum+1;
                        %     end
                        
                            ret = checkOverlap(P1,P2,Q1,Q2);
                            
                            if ret
                                L(i,1:6) = Li(idi,1:6);
                                L(i,7) = sqrt((Q1(1,1)- Q2(1,1))*(Q1(1,1) - Q2(1,1)) + (Q1(1,2) - Q2(1,2))*(Q1(1,2) - Q2(1,2)));
                                
                                %plot([CP(1,1) Mpi(idi,1)], [CP(1,2) Mpi(idi,2)],'-m','linewidth',1); hold on;
                                %plot([Q1(1,1) Q2(1,1)], [Q1(1,2) Q2(1,2)],'-g','linewidth',1); hold on;

                                notFound = 0;
                                F(i,1) = 1;
                            end                            
                        end
                        
                        %{
                        if (Fl == 1) % rotate and adjust line at k
                            F1(k,1) = 2;
                            rotrad = pi*rot/180;
                            cost = cos(rotrad); sint = sin(rotrad);
                            rotmat = [ cost sint; -sint cost];
                    
                            P1i = (rotmat*P1c')';
                            P2i = (rotmat*P2c')';
                            Pini = (rotmat*Pinc')';
                           
                            %indo = sid(k,1);
                            F(id(k,1),1) = 1; 
                            %V(id(k,1),1) = V(id(k,1),1)+1; 
                            P1r = [P1i(1,1) + CP(1,1) P1i(1,2) + CP(1,2)];
                            P2r = [P2i(1,1) + CP(1,1) P2i(1,2) + CP(1,2)];
                            Pinr = [Pini(1,1) + CP(1,1) Pini(1,2) + CP(1,2)];
                            L(id(k,1),1:2) = P1r;
                            L(id(k,1),3:4) = P2r;
                            L(id(k,1),8:9) = Pinr;
                        
                            L(id(k,1),5) = (P2r(1,2) - P1r(1,2))/(P2r(1,1) - P1r(1,1));
                            L(id(k,1),6) = P1r(1,2) -  L(id(k,1),5)*P1r(1,1);
                            
                            P1r = L(id(k,1),1:2);
                            P2r = L(id(k,1),3:4); 
                            %plot([P1r(1,1) P2r(1,1)], [P1r(1,2) P2r(1,2)],'-c','linewidth',1); hold on;
                            
                            %text(CP(1,1),CP(1,2),int2str(uint16(i))); hold on;
                            
                            adjust = 1;
                            D(id(k,1),1) = sqrt(cpsqr(now,1))/ResImage;
                        %else
                        %    V(id(k,1),1) = 1;
                        end %if Fl
                        %}
                    %end %if j
                    k = k+1;
                end%for k
                %now = now+1;
                %end%while
                
                %if adjust == 0
                %    V(id(k,1),1) = 1;
                %end
                %{
                if (adsum>0 && F(i,1)==0)
                    adsum = adsum+1;
                    %F1(j,1) = 2;
                    F(i,1) = 1; 
                    %V(id(j,1),1) = V(id(j,1),1)+1; 
                %else
                %    V(i,1) = V(i,1)+1; 
                %end%adsum
                total = total+adsum;
                %}
            %end %if F1
            %j = j+1;
            %count = count + adsum;
        %end%num
        %end%while count
        %always consider ref. line i as must selected curve
        %adsum = adsum+1;
        %F(i,1) = 1; 
        %total = total+adsum;
        %{
    elseif (V(i,1) == 0 && L(i,7) < minLength && F(i,1) == 0) %F(i,1) = 0: not visited yet; L(i,7) < minLength: this line cann't be a ref line to fix others with the same curveID; V(i,1) = 0: this line is by default selected
        %check the nearest ref. line (already selected)
        difX = Mp(:,1) - O*Mp(i,1);
        difY = Mp(:,2) - O*Mp(i,2);
        cp_sqr = difX.*difX + difY.*difY;
        
        [cp_sqr_sort id_sort] = sort(cp_sqr);
        refID = -1;
        for j=1:size(cp_sqr_sort,1)
            if (cp_sqr_sort(j,1) <= max_sqr_dist)
                if (L(id_sort(j,1),7) >= minLength)
                    refID = id_sort(j,1);
                    break;
                end
            else
                break;
            end
        end

        if (refID  > -1) % if a ref. found within max_buildingLength
            %add line i as non-ref line to be adjusted using the found ref
            %line
            %id = [i; id];
            %L1 = [L(i,:); L1];
            %F1 = [F(i,:); F1];
            %Mp1 = [Mp(i,:); Mp1];
            
            Q1 = L(refID,1:2);
                Q2 = L(refID,3:4);                
                 plot([Q1(1,1) Q2(1,1)], [Q1(1,2) Q2(1,2)],'-b','linewidth',2); hold on;
                 
                for k = 1:size(L1,1)                
                    if (refID ~= id(k,1) && F1(k,1) == 0)                     
                        CP = Mp1(k,:);
                        Q1c = [Q1(1,1) - CP(1,1) Q1(1,2) - CP(1,2)];
                        Q2c = [Q2(1,1) - CP(1,1) Q2(1,2) - CP(1,2)];
                        
                        m = (Q2c(1,1) - Q1c(1,1))/(Q2c(1,2) - Q1c(1,2));
                        theta = 180*atan(m)/pi;
                        mper = -1/m;
                        thetaper = 180*atan(mper)/pi;
                
                        P1 = L1(k,1:2);
                        P2 = L1(k,3:4); 
                        plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-r','linewidth',1); hold on;
                        Pin = L1(k,8:9);
                        P1c = [P1(1,1) - CP(1,1) P1(1,2) - CP(1,2)];
                        P2c = [P2(1,1) - CP(1,1) P2(1,2) - CP(1,2)];
                        Pinc = [Pin(1,1) - CP(1,1) Pin(1,2) - CP(1,2)];
                        
                        m1 = (P2c(1,1) - P1c(1,1))/(P2c(1,2) - P1c(1,2));
                        theta1 = 180*atan(m1)/pi;
                        thetadiff = abs(theta-theta1);
                        
                        Fl = 0;
                        if (thetadiff <= thetaThresh || abs(thetadiff-180) <= thetaThresh) % parallel to line at j of L1
                            rot = theta-theta1;
                            %L1(k,:) = rotateLine(L1(k,:), rot);
                            Fl = 1;
                            adsum = adsum+1;
                        else
                             thetadiff = abs(thetaper-theta1);
                             if (thetadiff <= thetaThresh || abs(thetadiff-180) <= thetaThresh) % perpendicular to line at j of L1
                                    rot = thetaper-theta1;
                                    %L1(k,:) = rotateLine(L1(k,:), rot);
                                    Fl = 1;
                                    adsum = adsum+1;
                             end
                        end
                        
                        if (Fl == 1) % rotate and adjust line at k
                            F1(k,1) = 2;
                            rotrad = pi*rot/180;
                            cost = cos(rotrad); sint = sin(rotrad);
                            rotmat = [ cost sint; -sint cost];
                    
                            P1i = (rotmat*P1c')';
                            P2i = (rotmat*P2c')';
                            Pini = (rotmat*Pinc')';
                           
                            %indo = sid(k,1);
                            F(id(k,1),1) = 1; 
                            %V(id(k,1),1) = V(id(k,1),1)+1; 
                            P1r = [P1i(1,1) + CP(1,1) P1i(1,2) + CP(1,2)];
                            P2r = [P2i(1,1) + CP(1,1) P2i(1,2) + CP(1,2)];
                            Pinr = [Pini(1,1) + CP(1,1) Pini(1,2) + CP(1,2)];
                            L(id(k,1),1:2) = P1r;
                            L(id(k,1),3:4) = P2r;
                            L(id(k,1),8:9) = Pinr;
                        
                            L(id(k,1),5) = (P2r(1,2) - P1r(1,2))/(P2r(1,1) - P1r(1,1));
                            L(id(k,1),6) = P1r(1,2) -  L(id(k,1),5)*P1r(1,1);
                            
                            %P1r = L(id(k,1),1:2);
                            %P2r = L(id(k,1),3:4); 
                            plot([P1r(1,1) P2r(1,1)], [P1r(1,2) P2r(1,2)],'-b','linewidth',1); hold on;
                        else
                            V(id(k,1),1) = 1;
                        end %if Fl
                    end %if j
                end%for k
            
        end
        %adsum = adsum+1;
        %F(i,1) = 1; 
        total = total+adsum;
        %}
    %end%if V
end%L
%hold off;

function ret = checkOverlap(P1,P2,Q1,Q2)
ret = 0;

%for P1P2
mp = (P2(1,2) - P1(1,2))/(P2(1,1) - P1(1,1));
cp = P2(1,2) - mp*P2(1,1);
mpper = -1/mp;
cpper1 = P1(1,2) - mpper*P1(1,1);
cpper2 = P2(1,2) - mpper*P2(1,1);

%for Q1Q2
mq = (Q2(1,2) - Q1(1,2))/(Q2(1,1) - Q1(1,1));
cq = Q2(1,2) - mq*Q2(1,1);
mqper = -1/mq;
cqper1 = Q1(1,2) - mqper*Q1(1,1);
cqper2 = Q2(1,2) - mqper*Q2(1,1);

%intersection points P1_|_Q1Q2
xi = (cq-cpper1)/(mpper-mq);
yi = mq*xi + cq;
%plot(xi,yi,'+r');
ci = yi - mqper*xi;

if sign(ci-cqper1) ~= sign(ci-cqper2)
    %plot(xi,yi,'+r');
    ret = 1;
end

if ret == 0
    %intersection points P2_|_Q1Q2
    xi = (cq-cpper2)/(mpper-mq);
    yi = mq*xi + cq;
    %plot(xi,yi,'+b');
    ci = yi - mqper*xi;    
    if sign(ci-cqper1) ~= sign(ci-cqper2)
        %plot(xi,yi,'+r');
        ret = 1;
    end
    
    if ret == 0;
        %intersection points Q1_|_P1P2
        xi = (cp-cqper1)/(mqper-mp);
        yi = mp*xi + cp;
        %plot(xi,yi,'+g');
        ci = yi - mpper*xi;
        if sign(ci-cpper1) ~= sign(ci-cpper2)
            %plot(xi,yi,'+r');
            ret = 1;
        end
        
        if ret == 0
            %intersection points Q2_|_P1P2
            xi = (cp-cqper2)/(mqper-mp);
            yi = mp*xi + cp;
            %plot(xi,yi,'+m');
            ci = yi - mpper*xi;
            if sign(ci-cpper1) ~= sign(ci-cpper2)
                %plot(xi,yi,'+r');
                ret = 1;
            end
        end
    end
end
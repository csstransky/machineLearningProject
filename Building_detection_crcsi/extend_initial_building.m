function [C1 P1] = extend_initial_building(A,R,mask,Rm,m,dAl,C1,P1,C2,ndvi,eMask,useEntropy,edgeThresh,ndviThresh,entropyThresh, ImgRect,FlagT)
%A,R,Img,nBin,m,dh,C1,P1,C2,col_hist,peaks,ALS,edgeThresh,htThresh,histThresh
                        m_left = -1/m;
                        c_left = P1(1,2) - m_left*P1(1,1);
                        %Y = m_left*C2(:,1) + c_left;
                        %Y_diff = Y-C2(:,2);
                        %len = sqrt((C1(1,1)-P1(1,1))*(C1(1,1)-P1(1,1)) + (C1(1,2)-P1(1,2))*(C1(1,2)-P1(1,2)));
                        %nPoints = dAc*dAl;
                        
                        %convert geographic to image
                            %for C1
                        %c1y = floor((C1(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - C1(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                        %c1x = floor((C1(1,2) - R(3,2) - c1y*R(2,2))/R(1,2));
                           %for P1
                        %p1y = floor((P1(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - P1(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                        %p1x = floor((P1(1,2) - R(3,2) - c1y*R(2,2))/R(1,2));
                        
                        %len = sqrt((c1x-p1x)*(c1x-p1x) + (c1y-p1y)*(c1y-p1y));
                        %count = 0;
                        %figure; mapshow(A,R); hold on;
                        while (1)                            
                            %if (count == 0)
                                Mp = [(C1(1,1)+P1(1,1))/2 (C1(1,2)+P1(1,2))/2];
                                ret1 = testInRectanglesAxisParallel(ImgRect(1,:), ImgRect(2,:), ImgRect(3,:), ImgRect(4,:), [C1]);
                                ret2 = testInRectanglesAxisParallel(ImgRect(1,:), ImgRect(2,:), ImgRect(3,:), ImgRect(4,:), [P1]);
                                retp = testInRectanglesAxisParallel(ImgRect(1,:), ImgRect(2,:), ImgRect(3,:), ImgRect(4,:), [Mp]);
                                if ((ret1 && retp) || (ret2 && retp))
                                    ret = 1;
                                else
                                    ret = 0;
                                end
                            %else
                            %    ret = testInRectanglesAxisParallel(ImgRect(1,:), ImgRect(2,:), ImgRect(3,:), ImgRect(4,:), [C1;P1]);
                            %end
                            if (ret)
                                [P1pa P2pa P3pa P4pa] = find_rectVertices(m_left,c_left,dAl,C1,P1);

                                %convert geographic to image
                                   %for P1
                                %p1pay = floor((P1pa(1,1)*R(1,2) - R(3,1)*R(1,2) + R(3,2)*R(1,1) - P1pa(1,2)*R(1,1))/(R(2,1)*R(1,2)-R(2,2)*R(1,1)));
                                %p1pax = floor((P1pa(1,2) - R(3,2) - c1y*R(2,2))/R(1,2));
                                
                                %wid = sqrt((c1x-p1pax)*(c1x-p1pax) + (c1y-p1pay)*(c1y-p1pay));
                                %nPoints = len*wid;
                                
                                %
                                
                                %plot([P1(1,1) C1(1,1)], [P1(1,2) C1(1,2)],'-y'); hold on; 
                                %plot([C2(1,1) C1(1,1)], [C2(1,2) C1(1,2)],'-b'); hold on; 
                                %plot(C2(1,1), C2(1,2),'*y'); plot(P1pa(1,1), P1pa(1,2),'+r'); plot(P2pa(1,1), P2pa(1,2),'ob');
                                %plot([P1pa(1,1) P2pa(1,1)], [P1pa(1,2) P2pa(1,2)],'-g'); hold on; 
                                %plot([P2pa(1,1) P4pa(1,1)], [P2pa(1,2) P4pa(1,2)],'-c'); hold on;
                                %plot([P3pa(1,1) P4pa(1,1)], [P3pa(1,2) P4pa(1,2)],'-y'); hold on;
                                %plot([P3pa(1,1) P1pa(1,1)], [P3pa(1,2) P1pa(1,2)],'-r'); hold on;
                                %

                                len1 = sqrt((C2(1,1)-P1pa(1,1))*(C2(1,1)-P1pa(1,1)) + (C2(1,2)-P1pa(1,2))*(C2(1,2)-P1pa(1,2)));
                                len2 = sqrt((C2(1,1)-P2pa(1,1))*(C2(1,1)-P2pa(1,1)) + (C2(1,2)-P2pa(1,2))*(C2(1,2)-P2pa(1,2)));

                                %Y1a = m_left*P1pa(:,1) + c_left;
                                %Y_diff1a = Y1a-P1pa(:,2);

                                %Y2a = m_left*P2pa(:,1) + c_left;
                                %Y_diff2a = Y2a-P2pa(:,2);

                                %if ((Y_diff <= 0 && Y_diff1a <= 0) || (Y_diff > 0 && Y_diff1a > 0)) % C2 and P1pa in same side                                
                                 %   ALS_left = find_ALS(A,R,ALS,P1,C1,P2pa,P4pa);
                                %else % C2 and P1pa in opposite sides                                
                                 %   ALS_left = find_ALS(A,R,ALS,P3pa,P1pa,C1,P1);
                                %end
                                
                                %if ((Y_diff <= 0 && Y_diff2a <= 0) || (Y_diff > 0 && Y_diff2a > 0)) % C2 and P2pa in same side                                
                                %mirP1 = mirror_point(m,c,P1);
                                
                                
                                C1m = obj2obs(C1,Rm);
                                P1m = obj2obs(P1,Rm);
                                
                                if (len1 > len2)
                                    
                                    P1pam = obj2obs(P1pa,Rm);
                                    %P2pam = obj2obs(P2pa,Rm);
                                    P3pam = obj2obs(P3pa,Rm);
                                    %P4pam = obj2obs(P4pa,Rm);
                                    
                                    %M_left = find_ALS(A,R,M,P3pa,P1pa,C1,P1);
                                    [totalaP numValaP] = findImagePointCounts(mask, P3pam, P1pam, C1m, P1m, 0);
                                    
                                    %mirP3pa = mirror_point(m,c,P3pa);
                                    %M_oth = find_ALS(A,R,M,mirP3pa,P1pa,C1,mirP1);                                
                                else % C2 and P2pa in opposite sides                                
                                    
                                    %P1pam = obj2obs(P1pa,Rm);
                                    P2pam = obj2obs(P2pa,Rm);
                                    %P3pam = obj2obs(P3pa,Rm);
                                    P4pam = obj2obs(P4pa,Rm);
                                    
                                    %M_left = find_ALS(A,R,M,P1,C1,P2pa,P4pa);
                                    [totalaP numValaP] = findImagePointCounts(mask, P1m,C1m,P2pam,P4pam, 0);
                                    
                                    %mirP4pa = mirror_point(m,c,P4pa);
                                    %M_oth = find_ALS(A,R,M,mirP1,C1,P2pa,mirP4pa);
                                end
                                %nPoints = max(size(M_left,1), size(M_oth,1));
                                %if (size(ALS_left,1) > 0)
                                     %ALSinfo = sum(ALS_left(:,3) >= htThresh)/size(ALS_left,1);
                                     %MInfo = sum((M_left(:,3) == 0))/nPoints;
                                     MInfo = numValaP/totalaP;%sum((M_left(:,3) == 0))/size(M_left,1);
                                     %MoInfo = sum((M_oth(:,3) == 0))/nPoints;
                                     %if (Flag == 1)
                                     %    F = MInfo >= edgeThresh & MoInfo < edgeThresh;
                                     %else
                                     %    F = MInfo >= edgeThresh;
                                     %end
                                     if (MInfo >= edgeThresh)
                                         C1i = obj2obs(C1,R);
                                         P1i = obj2obs(P1,R);
                                            
                                            if (len1 > len2)
                                                P1pai = obj2obs(P1pa,R);
                                                %P2pai = obj2obs(P2pa,R);
                                                P3pai = obj2obs(P3pa,R);
                                                %P4pai = obj2obs(P4pa,R);
                                                
                                                %N_left = find_ALS(A,R,N,P3pa,P1pa,C1,P1);
                                                m = findImagePointMean(ndvi, P3pai,P1pai,C1i,P1i);                                               
                                            else % C2 and P2pa in opposite sides                                
                                                %P1pai = obj2obs(P1pa,R);
                                                P2pai = obj2obs(P2pa,R);
                                                %P3pai = obj2obs(P3pa,R);
                                                P4pai = obj2obs(P4pa,R);
                                                
                                                %N_left = find_ALS(A,R,N,P1,C1,P2pa,P4pa);                                                                
                                                m = findImagePointMean(ndvi, P1i,C1i,P2pai,P4pai);
                                            end
                                            %m = mean(N_left(:,3));
                                            if (m <= ndviThresh && ~FlagT)
                                             %if ((Y_diff <= 0 && Y_diff1a <= 0) || (Y_diff > 0 && Y_diff1a > 0)) % C2 and P1pa in same side
                                             %   C1 = P2pa;
                                             %   P1 = P4pa;
                                            %else % C2 and P1pa in opposite sides
                                             %   C1 = P1pa;
                                                %P1 = P3pa;
                                             %end                                

                                              %if ((Y_diff <= 0 && Y_diff2a <= 0) || (Y_diff > 0 && Y_diff2a > 0)) % C2 and P1pa in same side
                                              if (len1 > len2)
                                                C1 = P1pa;
                                                P1 = P3pa;
                                              else % C2 and P1pa in opposite sides
                                                C1 = P2pa;
                                                P1 = P4pa;
                                              end        
                                            
                                            elseif (m > ndviThresh && FlagT && useEntropy)%if it is selected as a green roof edge
                                                if (len1 > len2)
                                                    %T_left = find_ALS(A,R,T,P3pa,P1pa,C1,P1);
                                                    [totaleP numValeP] = findImagePointCounts(eMask, P3pai,P1pai,C1i,P1i, 1);
                                                else % C2 and P2pa in opposite sides                                
                                                    %T_left = find_ALS(A,R,T,P1,C1,P2pa,P4pa);                                
                                                    [totaleP numValeP] = findImagePointCounts(eMask, P1i,C1i,P2pai,P4pai, 1);
                                                end
                                                m_entr = numValeP/totaleP;
                                                
                                                %id_edge = (T_left(:,3) == 1);
                                                %id_entr = (T_left(:,4) == 1);
                                                %m_edge = sum(id_edge)/size(T_left,1);
                                                %m_entr = sum(id_entr)/size(T_left,1);
                                                if (m_entr <= entropyThresh)% && m_edge > edgeShadowThresh)
                                                    if (len1 > len2)
                                                        C1 = P1pa;
                                                        P1 = P3pa;
                                                    else % C2 and P1pa in opposite sides
                                                        C1 = P2pa;
                                                        P1 = P4pa;
                                                    end 
                                                else
                                                    dAl = dAl/2;
                                                    %nPoints = (dAc*dAl)/(Sin*Sacross);
                                                    if (dAl < 0.4)
                                                        break;
                                                    end
                                                end
                                            else
                                                dAl = dAl/2;
                                                %nPoints = (dAc*dAl)/(Sin*Sacross);
                                                if (dAl < 0.4)
                                                    break;
                                                end
                                            end
                                     else
                                        dAl = dAl/2;
                                        %nPoints = (dAc*dAl)/(Sin*Sacross);
                                        if (dAl < 0.4)
                                            break;
                                        end
                                     end
                                %else
                                 %   break;
                                %end                            
                            else
                                break;
                            end
                        end
                        %hold off;
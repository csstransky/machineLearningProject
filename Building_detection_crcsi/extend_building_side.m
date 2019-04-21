function [C1 P1] = extend_building_side(A,R,m,dAl,C1,P1,C2,yiq3b,pMask,sMask,Rm,Mthresh,edgeThresh,HueThresh, SatThresh,IntThresh)
%function [C1 P1] = extend_segment_binary_image(A,R,m,dAc,dAl,C1,P1,C2,ALS,edgeThresh,Sin, Sacross)
%A,R,Img,nBin,m,dh,C1,P1,C2,col_hist,peaks,ALS,edgeThresh,htThresh,histThresh

dAlThresh = 0.1/abs(Rm(2,1));

                        m_left = -1/m;
                        c_left = P1(1,2) - m_left*P1(1,1);
                        %Y = m_left*C2(:,1) + c_left;
                        %Y_diff = Y-C2(:,2);
                        %dAc = sqrt((C1(1,1)-P1(1,1))*(C1(1,1)-P1(1,1)) + (C1(1,2)-P1(1,2))*(C1(1,2)-P1(1,2)));
                        %nPoints = (dAc*dAl)/(Sin*Sacross);
                        prev_MgInfo = 1;
                        MgFlag = 0;
                        while (1)                            
                            [P1pa P2pa P3pa P4pa] = find_rectVertices(m_left,c_left,dAl,C1,P1);
                            
                            %
                            %figure; mapshow(A,R); hold on;
                            %plot([P1(1,1) C1(1,1)], [P1(1,2) C1(1,2)],'-c'); hold on; 
                            %plot([C2(1,1) C1(1,1)], [C2(1,2) C1(1,2)],'-c'); hold on; 
                            %plot(C2(1,1), C2(1,2),'*c'); plot(P1pa(1,1), P1pa(1,2),'+r'); plot(P2pa(1,1), P2pa(1,2),'ob');
                            
                            %plot([P1pa(1,1) P2pa(1,1)], [P1pa(1,2) P2pa(1,2)],'-g'); hold on; 
                            %plot([P2pa(1,1) P4pa(1,1)], [P2pa(1,2) P4pa(1,2)],'-b'); hold on;
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
                            
                            C1m = obj2obs(C1,Rm);
                            P1m = obj2obs(P1,Rm);
                            
                            C1i = obj2obs(C1,R);
                            P1i = obj2obs(P1,R);
                                
                            if (len1 > len2)
                                P1pam = obj2obs(P1pa,Rm);
                                %P2pam = obj2obs(P2pa,Rm);
                                P3pam = obj2obs(P3pa,Rm);
                                %P4pam = obj2obs(P4pa,Rm);
                                
                                P1pai = obj2obs(P1pa,R);
                                %P2pai = obj2obs(P2pa,R);
                                P3pai = obj2obs(P3pa,R);
                                %P4pai = obj2obs(P4pa,R);
                                
                                %ALS_left = find_ALS(A,R,ALS,P3pa,P1pa,C1,P1);
                                ALS_left = findImagePoints(yiq3b, P3pai,P1pai,C1i,P1i);
                                
                                
                                %M_left = find_ALS(A,R,M,P3pa,P1pa,C1,P1);
                                [totalsP numValsP] = findImagePointCounts(sMask, P3pam, P1pam, C1m, P1m, 0);
                                
                                %Mg_left = find_ALS(A,R,Mg,P3pa,P1pa,C1,P1);
                                [totalpP numValpP] = findImagePointCounts(pMask, P3pam, P1pam, C1m, P1m, 0);
                                
                            else % C2 and P2pa in opposite sides                                
                                %P1pam = obj2obs(P1pa,Rm);
                                P2pam = obj2obs(P2pa,Rm);
                                %P3pam = obj2obs(P3pa,Rm);
                                P4pam = obj2obs(P4pa,Rm);
                                
                                %P1pai = obj2obs(P1pa,R);
                                P2pai = obj2obs(P2pa,R);
                                %P3pai = obj2obs(P3pa,R);
                                P4pai = obj2obs(P4pa,R);
                                
                                
                                %ALS_left = find_ALS(A,R,ALS,P1,C1,P2pa,P4pa);
                                ALS_left = findImagePoints(yiq3b, P1i,C1i,P2pai,P4pai);
                                
                                %M_left = find_ALS(A,R,M,P1,C1,P2pa,P4pa);
                                [totalsP numValsP] = findImagePointCounts(sMask, P1m,C1m,P2pam,P4pam, 0);
                                
                                %Mg_left = find_ALS(A,R,Mg,P1,C1,P2pa,P4pa);
                                [totalpP numValpP] = findImagePointCounts(pMask, P1m,C1m,P2pam,P4pam, 0);
                            end
                            nPoints = max([totalsP totalpP]);
                            if (size(ALS_left,1) > 0)
                                 %ALSinfo = sum(ALS_left(:,3) >= (meanHue - HueThresh) & ALS_left(:,3) <= (meanHue + HueThresh))/size(ALS_left,1);
                                 %ALSinfo = sum(ALS_left(:,3) >= meanHue)/size(ALS_left,1);
                                 
                                 HueBinary = logical(zeros(size(ALS_left,1),1));
                                 SatBinary = logical(zeros(size(ALS_left,1),1));
                                 IntBinary = logical(zeros(size(ALS_left,1),1));
                                    for j=1:size(HueThresh,1)
                                        BinH = (ALS_left(:,1) >= HueThresh(j,1)) & (ALS_left(:,1) < HueThresh(j,2)); 
                                        HueBinary = HueBinary | BinH;
                                    end
                                 
                                    for j=1:size(SatThresh,1)
                                        BinS = (ALS_left(:,2) >= SatThresh(j,1)) & (ALS_left(:,2) < SatThresh(j,2)); 
                                        SatBinary = SatBinary | BinS;
                                    end
                                    
                                    for j=1:size(IntThresh,1)
                                        BinI = (ALS_left(:,3) >= IntThresh(j,1)) & (ALS_left(:,3) < IntThresh(j,2)); 
                                        IntBinary = IntBinary | BinI;
                                    end
                                   
                                 %ALSinfo = sum(HueBinary,1)/nPoints;
                                 HueInfo = sum(HueBinary,1)/size(ALS_left,1);
                                 SatInfo = sum(SatBinary,1)/size(ALS_left,1);
                                 IntInfo = sum(IntBinary,1)/size(ALS_left,1);
                                 MInfo = numValsP/nPoints;%sum((M_left(:,3) == 0))/size(M_left,1);
                                 MgInfo = numValpP/nPoints;%sum((Mg_left(:,3) == 0))/size(Mg_left,1);
                                 
                                 if (MgFlag ==0 && MgInfo < 0.1)
                                     MgFlag = 1;
                                 end
                                 if (HueInfo >= edgeThresh && SatInfo >= edgeThresh && IntInfo >= edgeThresh && MInfo>= Mthresh)
                                     
                                     if (MgFlag == 0 || (MgFlag == 1 && prev_MgInfo >= MgInfo))
                                        prev_MgInfo = MgInfo;
                                          if (len1 > len2)
                                            C1 = P1pa;
                                            P1 = P3pa;
                                        else % C2 and P1pa in opposite sides
                                            C1 = P2pa;
                                            P1 = P4pa;
                                          end
                                     else
                                         break;
                                     end
                                 else
                                     dAl = dAl/2;
                                    %nPoints = (dAc*dAl)/(Sin*Sacross);
                                    if (dAl < dAlThresh)
                                        break;
                                    end
                                 end
                            else
                                break;
                            end
                        end
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
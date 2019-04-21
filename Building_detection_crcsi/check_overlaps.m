function ids = check_overlaps(A,R,Bground,max_buildingLength, P1, P2, P3, P4);

ids = [];

if (size(Bground,1)>0)
    for i = 1:size(Bground,1)
        p1 = Bground(i,1:2); 
        p3 = Bground(i,5:6);     
        cp(i,:) = [(p1(1,1) + p3(1,1))/2 (p1(1,2) + p3(1,2))/2];
    end

    CP = [(P1(1,1) + P3(1,1))/2 (P1(1,2) + P3(1,2))/2];

    O = ones(size(Bground,1),1);
    max_sqr_dist = max_buildingLength * max_buildingLength;
    difX = cp(:,1) - O*CP(1,1);
    difY = cp(:,2) - O*CP(1,2);
    cp_sqr = difX.*difX + difY.*difY;
        
    id = find(cp_sqr <= max_sqr_dist);
        %%%    
        %plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-g','LineWidth', 1); hold on;
        %plot([P2(1,1) P3(1,1)], [P2(1,2) P3(1,2)],'-g','LineWidth', 1); hold on;
        %plot([P3(1,1) P4(1,1)], [P3(1,2) P4(1,2)],'-g','LineWidth', 1); hold on;
        %plot([P4(1,1) P1(1,1)], [P4(1,2) P1(1,2)],'-g','LineWidth', 1); hold on;
        %%%
        for k = 1:size(id,1)
            j = id(k,1);
            %if (i ~= j && Flag(j,1) == 1)
                Pj1 = Bground(j,1:2); Pj2 = Bground(j,3:4); Pj3 = Bground(j,5:6); Pj4 = Bground(j,7:8);
                
                mp12_1 = [(Pj1(1,1) + Pj2(1,1))/2 (Pj1(1,2) + Pj2(1,2))/2];
                mp12_2 = [(Pj1(1,1) + mp12_1(1,1))/2 (Pj1(1,2) + mp12_1(1,2))/2];
                mp12_3 = [(mp12_1(1,1) + Pj2(1,1))/2 (mp12_1(1,2) + Pj2(1,2))/2];
                %plot(mp12_1(1,1),mp12_1(1,2),'+g'); plot(mp12_2(1,1),mp12_2(1,2),'+b'); plot(mp12_3(1,1),mp12_3(1,2),'+r');
                
                mp23_1 = [(Pj3(1,1) + Pj2(1,1))/2 (Pj3(1,2) + Pj2(1,2))/2];
                mp23_2 = [(Pj3(1,1) + mp23_1(1,1))/2 (Pj3(1,2) + mp23_1(1,2))/2];
                mp23_3 = [(mp23_1(1,1) + Pj2(1,1))/2 (mp23_1(1,2) + Pj2(1,2))/2];
                %plot(mp23_1(1,1),mp23_1(1,2),'+g'); plot(mp23_2(1,1),mp23_2(1,2),'+b'); plot(mp23_3(1,1),mp23_3(1,2),'+r');
                
                mp34_1 = [(Pj3(1,1) + Pj4(1,1))/2 (Pj3(1,2) + Pj4(1,2))/2];
                mp34_2 = [(Pj3(1,1) + mp34_1(1,1))/2 (Pj3(1,2) + mp34_1(1,2))/2];
                mp34_3 = [(mp34_1(1,1) + Pj4(1,1))/2 (mp34_1(1,2) + Pj4(1,2))/2];
                % plot(mp34_1(1,1),mp34_1(1,2),'+g'); plot(mp34_2(1,1),mp34_2(1,2),'+b'); plot(mp34_3(1,1),mp34_3(1,2),'+r');
                 
                mp41_1 = [(Pj1(1,1) + Pj4(1,1))/2 (Pj1(1,2) + Pj4(1,2))/2];
                mp41_2 = [(Pj1(1,1) + mp41_1(1,1))/2 (Pj1(1,2) + mp41_1(1,2))/2];
                mp41_3 = [(mp41_1(1,1) + Pj4(1,1))/2 (mp41_1(1,2) + Pj4(1,2))/2];
                % plot(mp41_1(1,1),mp41_1(1,2),'+g'); plot(mp41_2(1,1),mp41_2(1,2),'+b'); plot(mp41_3(1,1),mp41_3(1,2),'+r');
                 
                %ret = testInRectangles(P1, P2, P3, P4, [Pj1;Pj2;Pj3;Pj4;mp12_1;mp12_2;mp12_3;mp23_1;mp23_2;mp23_3;mp34_1;mp34_2;mp34_3;mp41_1;mp41_2;mp41_3]);
                %ret = testInRectangles_number(P1, P2, P3, P4, [Pj1;Pj2;Pj3;Pj4]);
                ret = testInRectangles_number(Pj1,Pj2,Pj3,Pj4, [P1; P2; P3; P4]);
                
                 %%%    
                %plot([Pj1(1,1) Pj2(1,1)], [Pj1(1,2) Pj2(1,2)],'-b','LineWidth', 1); hold on;
                %plot([Pj2(1,1) Pj3(1,1)], [Pj2(1,2) Pj3(1,2)],'-b','LineWidth', 1); hold on;
                %plot([Pj3(1,1) Pj4(1,1)], [Pj3(1,2) Pj4(1,2)],'-b','LineWidth', 1); hold on;
                %plot([Pj4(1,1) Pj1(1,1)], [Pj4(1,2) Pj1(1,2)],'-b','LineWidth', 1); hold on;
                %%%

                if (ret > 0)                    
                    ids = [ids; [j ret]]; 
                end            
        end
end
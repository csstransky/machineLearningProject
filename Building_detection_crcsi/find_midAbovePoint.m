function P = find_midAbovePoint(midP,m,c,matP);
% parallel line going through matP
m1 = m;
c1 = matP(1,2) - m1*matP(1,1);

%perpendicular line through midP
m2 = -1/m;
c2 = midP(1,2) - m2*midP(1,1);

P(1,1) = (c2-c1)/(m1-m2);
P(1,2) = m1*P(1,1)+c1;
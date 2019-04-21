function mP = mirror_point(m,c,P);
m1 = -1/m;
c1 = P(1,2)-m1*P(1,1);

%find intersection point
intP(1,1) = (c-c1)/(m1-m);
intP(1,2) = m*intP(1,1)+c;

%find mirror point
mP(1,1) = 2*intP(1,1) - P(1,1);
mP(1,2) = 2*intP(1,2) - P(1,2);
function Pout = obj2obs(P,Rd)

y1 = floor((P(1,1)*Rd(1,2) - Rd(3,1)*Rd(1,2) + Rd(3,2)*Rd(1,1) - P(1,2)*Rd(1,1))/(Rd(2,1)*Rd(1,2)-Rd(2,2)*Rd(1,1)));
x1 = floor((P(1,2) - Rd(3,2) - y1*Rd(2,2))/Rd(1,2));

idy = y1 <= 0; y1(idy,:) = 1;
idx = x1 <= 0; x1(idx,:) = 1;

Pout = [x1 y1];
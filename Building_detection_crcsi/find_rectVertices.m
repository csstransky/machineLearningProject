function [P1 P2 P3 P4] = find_rectVertices(m,c,d,C1,C2);

%     P1-----------P3
%     |                   |d
%    C1----------C2
 %    |                   |d
%     P2----------P4

co1 = C1(1,2) + C1(1,1)/m;  % y-intersection of first perpendicular line
co2 = C2(1,2) + C2(1,1)/m;  % y-intersection of second perpendicular line

e = 1 + (1/(m*m));
mde = -2*e;
sd = d*d;

% coeffiecints of first quadratic equeation to solve
A = e;
B = mde*C1(1,1);
C = C1(1,1)*C1(1,1)*e - sd;

dA = 2*A;
Det = sqrt(B*B - 4*A*C);

P1(1,1) = (-B + Det)/dA;
P1(1,2) = -P1(1,1)/m + co1;

P2(1,1) = (-B - Det)/dA;
P2(1,2) = -P2(1,1)/m + co1;

% coeffiecints of first quadratic equeation to solve
%A = e;
B = mde*C2(1,1);
C = C2(1,1)*C2(1,1)*e - sd;

%dA = A*A;
Det = sqrt(B*B - 4*A*C);

P3(1,1) = (-B + Det)/dA;
P3(1,2) = -P3(1,1)/m + co2;

P4(1,1) = (-B - Det)/dA;
P4(1,2) = -P4(1,1)/m + co2;

%figure;
%plot([P1(1,1) P2(1,1)], [P1(1,2) P2(1,2)],'-g'); hold on;
%plot([P2(1,1) P4(1,1)], [P2(1,2) P4(1,2)],'-b'); hold on;
%plot([P4(1,1) P3(1,1)], [P4(1,2) P3(1,2)],'-y'); hold on;
%plot([P3(1,1) P1(1,1)], [P3(1,2) P1(1,2)],'-r'); hold on;
%plot([C1(1,1) C2(1,1)], [C1(1,2) C2(1,2)],'-c'); hold on;
%here = 1;
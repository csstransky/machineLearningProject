function yes = isInsideRectangle(P1, P2, P3, P4, P)

	% calculate (twice, multiplied by 2) areas of clockwise triangles with P
	% for clockwise areas (of all 4) triangles (P1P2P, P2P3P, P3P4P, P4P1P) should be -ve, if P is inside the quadrangle
	% for anti-clockwise all areas should be +ve, if P is inside the quadrangle
	a1 = P2(1,1)*P(1,2) - P2(1,2)*P(1,1) - P1(1,1)*P(1,2) + P1(1,2)*P(1,1) + P1(1,1)*P2(1,2) - P1(1,2)*P2(1,1);
	a2 = P3(1,1)*P(1,2) - P3(1,2)*P(1,1) - P2(1,1)*P(1,2) + P2(1,2)*P(1,1) + P2(1,1)*P3(1,2) - P2(1,2)*P3(1,1);
	a3 = P4(1,1)*P(1,2) - P4(1,2)*P(1,1) - P3(1,1)*P(1,2) + P3(1,2)*P(1,1) + P3(1,1)*P4(1,2) - P3(1,2)*P4(1,1);
	a4 = P1(1,1)*P(1,2) - P1(1,2)*P(1,1) - P4(1,1)*P(1,2) + P4(1,2)*P(1,1) + P4(1,1)*P1(1,2) - P4(1,2)*P1(1,1);
	
    yes = (a1<0 && a2<0 && a3<0 && a4<0) || (a1>0 && a2>0 && a3>0 && a4>0);
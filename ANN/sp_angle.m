function Output = sp_angle(c_A, c_B)
% find super patch angle base on given centers barycenters
x = c_B(1) - c_A(1);
y = c_B(2) - c_A(2);

angle = atan2(y, x);

if (c_B(2) < c_A(2))
    Output = -angle;
else
    Output = 2*pi - angle;
end

end


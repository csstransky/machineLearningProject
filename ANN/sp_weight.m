function Output = sp_weight(c_Ai, c_Ai_, c_Bj, c_Bj_, m, n, nA, radius)

ws_Ai_ = exp(-((c_Ai(2)-c_Ai_(2))^2+(c_Ai(1)-c_Ai_(1))^2)/(2*radius^2));
ws_Bj_ = exp(-((c_Bj(2)-c_Bj_(2))^2+(c_Bj(1)-c_Bj_(1))^2)/(2*radius^2));

Output = exp(-(c_Bj_ - c_Ai_+ c_Ai - c_Bj)*(c_Bj_-c_Ai_+c_Ai-c_Bj)'/(0.5 ...
    *sqrt(m*n/nA))) * ws_Ai_ * ws_Bj_;
end


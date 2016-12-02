function [x_params, analytic, sampled, pem_direct, pem_iter, hong, hongs, zhao, points, weights] = pem(yfun, x, x_2, pdf_x, type, a, b)
        v = skewness(x);
        e1 = v/2 + sqrt(1+(v/2)^2);
        e2 = e1 - v;
        p1 = e2/(e1+e2);
        p2 = 1-p1;
        x1 = mean(x) + std(x)*e1;
        x2 = mean(x) - std(x)*e2;
    
        pem_mean = p1*yfun(x1) + p2*yfun(x2);
        pem_sig = sqrt(p1*p2)*abs(yfun(x1) - yfun(x2));
        pem_sig_iter = sqrt(p1*(yfun(x1) - pem_mean)^2 + p2*(yfun(x2) - pem_mean)^2);
        pem_v_y = (1/pem_sig_iter)*(p2 - p1)*(yfun(x1) - yfun(x2));
        %Modify PEM here
        %Idea = re-weight P1 and P2 according to distance from t he mean of
        %Y, which is now know
        pem_v_y_iter = (((yfun(x1) - pem_mean)/pem_sig)^3 + ((yfun(x2) - pem_mean)/pem_sig)^3)/2;
        
        mean_actual = mean(yfun(x));
        sig_actual = std(yfun(x));
        skew_actual = skewness(yfun(x));
        
        points = [x1, x2];
        weights = [p1, p2];
        
        %hongs pem for k=2:
        e1 = v/2 + sqrt(1 + (v/2)^2);
        e2 = v/2 - sqrt(1 + (v/2)^2);
        p1 = -e2/(e1-e2);
        p2 = e1/(e1-e2);
        x1 = mean(x) + e1*std(x);
        x2 = mean(x) + e2*std(x);
        hong_mean = p1*yfun(x1) + p2*yfun(x2);
        hong_std_raw = p1*(yfun(x1))^2 + p2*(yfun(x2))^2;
        hong_skew_raw = (p1*(yfun(x1))^3 + p2*(yfun(x2))^3);
        hong_std = sqrt(hong_std_raw - hong_mean^2);
        hong_third_mom = hong_skew_raw -3*hong_mean*hong_std_raw + 2*hong_mean^3;
        hong_skew = hong_third_mom/hong_std^3;
        points = [points, x1, x2];
        weights = [weights, p1, p2];
        
        %hongs pem k = 3
        e1 = v/2 + sqrt(kurtosis(x) - (3/4)*v^2);
        e2 = v/2 - sqrt(kurtosis(x) - (3/4)*v^2);
        e3 = 0;
        w1 = 1/(e1*(e1 - e2));
        w2 = -1/(e2*(e1-e2));
        w3 = 1 - (1/(kurtosis(x) - v^2));
        
        x1 = mean(x) + e1*std(x);
        x2 = mean(x) + e2*std(x);
        x3 = mean(x) + e3*std(x);
        points = [points, x1,x2,x3];
        weights = [weights, w1,w2,w3];
        
        hongs_mean = w1*yfun(x1) + w2*yfun(x2) + w3*yfun(x3);
        hongs_std_raw = w1*(yfun(x1))^2 + w2*(yfun(x2))^2 + w3*yfun(x3)^2;
        hongs_skew_raw = (w1*(yfun(x1))^3 + w2*(yfun(x2))^3 + w3*(yfun(x3))^3);
        hongs_std = sqrt(hongs_std_raw - hongs_mean^2);
        hongs_third_mom = hongs_skew_raw -3*hongs_mean*hongs_std_raw + 2*hongs_mean^3;
        hongs_skew = hongs_third_mom/hongs_std^3;
        if type == 'U'
            analytic_mean = a*sqrt(pi/2);
            analytic_std = sqrt((4-pi)*a^2/2);
            analytic_skewness = (2*sqrt(pi)*(pi-3))/(4-pi)^1.5;
        end
        if type == 'R'
            analytic_mean = 2;
            analytic_std = 2;
            analytic_skewness = 2;
        end
        if type == 'B'
            analytic_mean = a/(a+b) + 2;
            analytic_std = sqrt((a*b)/((a+b+1)*(a+b)^2));
            analytic_skewness = (2*(b - a)*sqrt(a+b+1))/((a+b+2)*sqrt(a*b));
        end
        
        %zhaos pem
        cdf_x = cumtrapz(x_2, pdf_x);
        non_mono_indices = find(cdf_x == 0);
        non_mono_indices = non_mono_indices(1:end - 1);
        non_mono_indices = [non_mono_indices, find(diff(cdf_x) == 0)];
        cdf_x(non_mono_indices) = [];
        x_2(non_mono_indices) = [];

        p0 = 16/35;
        u1pos = 1.1544054;
        p1 = 0.2401233;
        u1neg = -u1pos;
        u2pos = 2.3667594;
        u2neg = -u2pos;
        p2 = 3.07571*10^-2;
        u3pos = 3.7504397;
        u3neg = -u3pos;
        p3 = 5.48269*10^-4;

        x0 = cdf('Normal', 0, 0, 1);
        x0 = interp1(cdf_x,x_2, x0)

        x1pos = cdf('Normal', u1pos, 0, 1);
        x1pos = interp1(cdf_x , x_2 , x1pos)

        x1neg = cdf('Normal', u1neg, 0, 1);
        x1neg = interp1(cdf_x , x_2 , x1neg)

        x2pos = cdf('Normal', u2pos, 0, 1);
        x2pos = interp1(cdf_x , x_2 , x2pos)

        x2neg = cdf('Normal', u2neg, 0, 1);
        x2neg = interp1(cdf_x , x_2 , x2neg)

        x3pos = cdf('Normal', u3pos, 0, 1);
        x3pos = interp1(cdf_x , x_2 , x3pos)

        x3neg = cdf('Normal', u3neg, 0, 1);
        x3neg = interp1(cdf_x , x_2 , x3neg)

        zhao_mean = p0*yfun(x0) + p1*yfun(x1pos) + p1*yfun(x1neg) + p2*yfun(x2pos) + p2*yfun(x2neg) + p3*yfun(x3neg) + p3*yfun(x3pos);

        zhao_std = sqrt(p0*(yfun(x0) - zhao_mean)^2 + p1*(yfun(x1pos) - zhao_mean)^2 + p1*(yfun(x1neg) - zhao_mean)^2 + p2*(yfun(x2pos) - zhao_mean)^2 + p2*(yfun(x2neg) - zhao_mean)^2 + p3*(yfun(x3neg) - zhao_mean)^2 + p3*(yfun(x3pos)- zhao_mean)^2);

        zhao_skew = (1/zhao_std^3)*(p0*(yfun(x0) - zhao_mean)^3 + p1*(yfun(x1pos) - zhao_mean)^3 + p1*(yfun(x1neg) - zhao_mean)^3 + p2*(yfun(x2pos) - zhao_mean)^3 + p2*(yfun(x2neg) - zhao_mean)^3 + p3*(yfun(x3neg) - zhao_mean)^3 + p3*(yfun(x3pos)- zhao_mean)^3);
        
        points = [points, x0, x1pos, x1neg, x2pos, x2neg, x3pos, x3neg];
        weights = [weights, p0, p1, p2, p3];
        %zhaos pem ends 
        
        pem_direct = [pem_mean; pem_sig; pem_v_y];
        pem_iter = [pem_mean; pem_sig_iter; pem_v_y_iter];

        sampled = [mean_actual; sig_actual; skew_actual];

        analytic = [analytic_mean; analytic_std; analytic_skewness];

        x_params = [mean(x); std(x); skewness(x)];
        
        hong = [hong_mean; hong_std; hong_skew];
        
        hongs = [hongs_mean; hongs_std; hongs_skew];
        
        zhao = [zhao_mean; zhao_std; zhao_skew];
        
        
       
end
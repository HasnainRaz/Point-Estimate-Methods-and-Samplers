clc;
clear all;
close all;
rng(9)

n = 25000;
n2 = 20000;

s_tech = 'lhs';
dist = 'norm';

if strcmp(dist, 'norm')
    x0 = 200:0.1:300;
    x1 = 210:0.1:310;
    x2 = 30:0.1:180;
    x3 = 50:0.1:200;

    pd1 = normpdf(x0, 225, 3);
    pd2 = normpdf(x1, 265, 3);
    pd3 = normpdf(x2, 65, 3);
    pd4 = normpdf(x3, 85, 3);
    
elseif strcmp(dist, 'beta')
    
    x_all = betarnd(2,5, 1, 10000);
    x0 = (235 - 215)*x_all + 215;
    x1 = (275 - 255)*x_all + 255;
    x2 = (75 - 55)*x_all + 55;
    x3 = (95 - 75)*x_all + 75;
    
    [pd1, ~] = ksdensity(x0, x0);
    [pd2, ~] = ksdensity(x1, x1);
    [pd3, ~] = ksdensity(x2, x2);
    [pd4, ~] = ksdensity(x3, x3);
end

if strcmp(s_tech, 'mh')
    [pr1, ~, ~] = samplers(s_tech, x0, pd1, n);
    [pr2, ~, ~] = samplers(s_tech, x0, pd1, n);
    [pr3, ~, ~] = samplers(s_tech, x1, pd2, n);
    
    [pi1, ~, ~] = samplers(s_tech, x2, pd3, n);
    [pi2, ~, ~] = samplers(s_tech, x2, pd3, n);
    [pi3, ~, ~] = samplers(s_tech, x3, pd4, n);   
else
    [pr1, ~, ~] = samplers(s_tech, x0, pd1, n);
    [pr2, ~, ~] = samplers(s_tech, x0, pd1, n);
    [pr3, ~, ~] = samplers(s_tech, x1, pd2, n);
    
    [pi1, ~, ~] = samplers(s_tech, x2, pd3, n);
    [pi2, ~, ~] = samplers(s_tech, x2, pd3, n);
    [pi3, ~, ~] = samplers(s_tech, x3, pd4, n);
end

pr1 = pr1(randperm(length(pr1)));
pr2 = pr2(randperm(length(pr2)));
pr3 = pr3(randperm(length(pr3)));
pi1 = pi1(randperm(length(pi1)));
pi2 = pi2(randperm(length(pi2)));
pi3 = pi3(randperm(length(pi3)));

mopt = mpoption('verbose', 0, 'out.all', 0);
result = runpf('case5', mopt);

tic
for i=1:1:n2
    result.bus(2, 3) = pr1(i);
    result.bus(3, 3) = pr2(i);
    result.bus(4, 3) = pr3(i);
    
    result.bus(2, 4) = pi1(i);
    result.bus(3, 4) = pi2(i);
    result.bus(4, 4) = pi3(i);
    
    result = runpf(result, mopt);
    
    gen_real(1:4, i) = result.gen(:, 2);
    gen_reac(1:4, i) = result.gen(:, 3);
end
toc

%Change this to only consider real and im power variables, hence, only one
%or 2 dimensional cases

figure
histogram(gen_real, 'normalization', 'pdf')
title('Total real power generation')

tic
%point estimate here:
result = runpf('case5', mopt);
if strcmp(dist, 'norm')
    x1 = normrnd(225, 3, 1, 150);
    x2 = normrnd(265, 3, 1, 150);
    x3 = normrnd(65, 3, 1, 150);
    x4 = normrnd(85, 3, 1, 150);
elseif strcmp(dist, 'beta')
    x_all = betarnd(2,5, 1, 10000);
    x1 = (235 - 215)*x_all + 215;
    x2 = (275 - 255)*x_all + 255;
    x3 = (75 - 55)*x_all + 55;
    x4 = (95 - 75)*x_all + 75;
end

m1 = 225;
m2 = 225;
m3 = 265;
mi1 = 65;
mi2 = 65;
mi3 = 85;

[p1, w1] = pem_mini('hong3', x1);
[p2, w2] = pem_mini('hong3', x2);
[p3, w3] = pem_mini('hong3', x3);
[p4, w4] = pem_mini('hong3', x4);

result.bus(2, 3) = p1(1);
result.bus(3, 3) = m2;
result.bus(4, 3) = m3;

result.bus(2, 4) = mi1;
result.bus(3, 4) = mi2;
result.bus(4, 4) = mi3;

result = runpf(result, mopt);
pem_real_gen(1:4, 1) = result.gen(:, 2);

result.bus(2, 3) = p1(2);
result.bus(3, 3) = m2;
result.bus(4, 3) = m3;

result.bus(2, 4) = mi1;
result.bus(3, 4) = mi2;
result.bus(4, 4) = mi3;

result = runpf(result, mopt);
pem_real_gen(1:4, 2) = result.gen(:, 2);

result.bus(2, 3) = p1(3);
result.bus(3, 3) = m2;
result.bus(4, 3) = m3;

result.bus(2, 4) = mi1;
result.bus(3, 4) = mi2;
result.bus(4, 4) = mi3;

result = runpf(result, mopt);
pem_real_gen(1:4, 3) = result.gen(:, 2);

hong_mean = w1(1).*pem_real_gen(: , 1) + w1(2).*pem_real_gen(:, 2) + w1(3).*pem_real_gen(:, 3);

hong_std_raw = w1(1).*pem_real_gen(:, 1).^2 + w1(2).*pem_real_gen(:, 2).^2 + w1(3).*pem_real_gen(:, 3).^2;

hong_std = sqrt(hong_std_raw - hong_mean.^2);

hong_skew_raw = w1(1)*pem_real_gen(:, 1).^3 + w1(2)*pem_real_gen(:, 2).^3 + w1(3)*pem_real_gen(:, 3);

hong_third_mom = (hong_skew_raw - 3*hong_mean.*hong_std_raw + 2*hong_mean.^3)./hong_std.^3;

disp('True Values: ')
true_mean = [mean(gen_real(1, :)); mean(gen_real(2, :)); mean(gen_real(3, :)); mean(gen_real(4, :))];

true_std = [std(gen_real(1, :)); std(gen_real(2, :)); std(gen_real(3, :)); std(gen_real(4, :))];

true_skew = [skewness(gen_real(1, :)); skewness(gen_real(2, :)); skewness(gen_real(3, :)); skewness(gen_real(4, :))];

table(true_mean, hong_mean, true_std, hong_std, true_skew, hong_third_mom, 'RowNames', {'Gen_1', 'Gen_2', 'Gen_3', 'Gen_4'}, 'VariableNames', {'True_Mean', 'Hong_Mean', 'True_Std', 'Hong_Std', 'True_Skew', 'Hong_Skew'})


% result.bus(2, 3) = m1;
% result.bus(3, 3) = p1(1);
% result.bus(4, 3) = m3;
% 
% result.bus(2, 4) = mi1;
% result.bus(3, 4) = mi2;
% result.bus(4, 4) = mi3;
% 
% result = runpf(result, mopt);
% pem_real_gen(3) = sum(result.gen(:, 2));
% second(3) = sum(result.gen(:, 3));
% 
% result.bus(2, 3) = m1;
% result.bus(3, 3) = p1(2);
% result.bus(4, 3) = m3;
% 
% result.bus(2, 4) = mi1;
% result.bus(3, 4) = mi2;
% result.bus(4, 4) = mi3;
% 
% result = runpf(result, mopt);
% pem_real_gen(4) = sum(result.gen(:, 2));
% second(4) = sum(result.gen(:, 3));
% 
% result.bus(2, 3) = m1;
% result.bus(3, 3) = m2;
% result.bus(4, 3) = p2(1);
% 
% result.bus(2, 4) = mi1;
% result.bus(3, 4) = mi2;
% result.bus(4, 4) = mi3;
% 
% result = runpf(result, mopt);
% pem_real_gen(5) = sum(result.gen(:, 2));
% second(5) = sum(result.gen(:, 3));
% 
% result.bus(2, 3) = m1;
% result.bus(3, 3) = m2;
% result.bus(4, 3) = p2(2);
% 
% result.bus(2, 4) = mi1;
% result.bus(3, 4) = mi2;
% result.bus(4, 4) = mi3;
% 
% result = runpf(result, mopt);
% pem_real_gen(6) = sum(result.gen(:, 2));
% second(6) = sum(result.gen(:, 3));
% 
% mr1 = pem_real_gen(1)*w1(1) + pem_real_gen(2)*w1(2);
% mr2 = pem_real_gen(3)*w1(1) + pem_real_gen(4)*w1(2);
% mr3 = pem_real_gen(5)*w2(1) + pem_real_gen(6)*w2(2);
% means = [mr1, mr2, mr3];
% 
% hong_std_raw(1) = w1(1)*pem_real_gen(1)^2 + w1(2)*pem_real_gen(2)^2;
% hong_std_raw(2) = w1(1)*pem_real_gen(3)^2 + w1(2)*pem_real_gen(4)^2;
% hong_std_raw(3) = w2(1)*pem_real_gen(5)^2 + w2(2)*pem_real_gen(6)^2;
% 
% hong_skew_raw(1) = (w1(1)*pem_real_gen(1)^3 + w1(2)*pem_real_gen(2)^3);
% hong_skew_raw(2) = (w1(1)*pem_real_gen(3)^3 + w1(2)*pem_real_gen(4)^3);
% hong_skew_raw(3) = (w2(1)*pem_real_gen(5)^3 + w2(2)*pem_real_gen(6)^3);
% 
% hong_std(1) = sqrt(hong_std_raw(1) - means(1)^2);
% hong_std(2) = sqrt(hong_std_raw(2) - means(2)^2);
% hong_std(3) = sqrt(hong_std_raw(3) - means(3)^2);
% 
% hong_third_mom(1) = hong_skew_raw(1) -3*means(1)*hong_std_raw(1) + 2*means(1)^3;
% hong_third_mom(2) = hong_skew_raw(2) -3*means(2)*hong_std_raw(2) + 2*means(2)^3;
% hong_third_mom(3) = hong_skew_raw(3) -3*means(3)*hong_std_raw(3) + 2*means(3)^3;
% 
% hong_skew(1) = hong_third_mom(1)/hong_std(1)^3;
% hong_skew(2) = hong_third_mom(2)/hong_std(2)^3;
% hong_skew(3) = hong_third_mom(3)/hong_std(3)^3;
% 
% result.bus(2, 3) = m1;
% result.bus(3, 3) = m2;
% result.bus(4, 3) = m3;
% 
% result.bus(2, 4) = mi1;
% result.bus(3, 4) = mi2;
% result.bus(4, 4) = mi3;
% 
% result = runpf(result, mopt);
% mean_mr = sum(result.gen(:, 2));
% mean_mi = sum(result.gen(:, 3));
% toc
% 
% disp('PEM Std Real Power Generation')
% mean(hong_std)
% 
% disp('PEM Skew Real Power Generation')
% mean(hong_skew)
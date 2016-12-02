clc;
clear all;
close all;

n = 100000;
rows = {'mean', 'standard dev', 'skewness'};
cols = {'X', 'Y_analytic', 'Y_sampled', 'Rosen_pem_direct', 'Rosen_pem_modded', 'Hong2', 'Hong3', 'Zhao'};

a = 1;
x = rand(1,n);
x_2 = -4:0.0001:8;
pdf_x = unifpdf(x_2, 0,1);
yfun = @(x)a*sqrt(-2*log(1-x));
[x_params, analytic, sampled, pem_direct, pem_iter, hong, hongs, zhao, points, weights] = pem(yfun, x, x_2, pdf_x, 'U', a);
disp(['X = Uniform, Y = Rayleigh(', num2str(a), ')'])
weights(1:4)
points(1:4)
table(x_params, analytic, sampled, pem_direct, pem_iter, hong, hongs, zhao, 'RowNames', rows, 'VariableNames', cols)
figure()
hold on
plot(-1:0.001:2, unifpdf(-1:0.001:2, 0,1))
p1 = plot([points(1), points(1)], [0, weights(1)], 'g');
plot([points(2), points(2)], [0, weights(2)], 'g')
p2 = plot([points(5), points(5)], [0, weights(5)], 'b');
plot([points(6), points(6)], [0, weights(6)], 'b')
plot([points(7), points(7)], [0, weights(7)], 'b')

p3 = plot([points(8), points(8)], [0, weights(8)], 'r');
plot([points(9), points(9)], [0, weights(9)], 'r')
plot([points(10), points(10)], [0, weights(9)], 'r')
plot([points(11), points(11)], [0, weights(10)], 'r')
plot([points(12), points(12)], [0, weights(10)], 'r')
plot([points(13), points(13)], [0, weights(11)], 'r')
plot([points(14), points(14)], [0, weights(11)], 'r')
legend([p1,p2,p3], 'Rosen/Hong 2 point', 'Hong 3 point', 'Zhao 7 point')
hold off
title('X = Uniform')
ylim([0, 1.2])
y = raylpdf(0:0.1:4, a);
figure()
plot(0:0.1:4,y)
title(['Y = Rayleigh(',num2str(a),')'])



x = raylrnd(1, 1, n);
yfun = @(x)x.*x;
pdf_x = raylpdf(x_2,1);
[x_params, analytic, sampled, pem_direct, pem_iter, hong, hongs, zhao, points, weights] = pem(yfun, x, x_2, pdf_x, 'R');
disp('X = Rayleigh(1), Y = Chi-squared(2)')
table(x_params, analytic, sampled, pem_direct, pem_iter, hong, hongs, zhao, 'RowNames', rows, 'VariableNames', cols)
weights(1:4)
points(1:4)
figure()
hold on
plot(0:0.1:5, raylpdf(0:0.1:5, 1))
p1 = plot([points(1), points(1)], [0, weights(1)], 'g');
plot([points(2), points(2)], [0, weights(2)], 'g')

p2 = plot([points(5), points(5)], [0, weights(5)], 'b');
plot([points(6), points(6)], [0, weights(6)], 'b')
plot([points(7), points(7)], [0, weights(7)], 'b')

p3 = plot([points(8), points(8)], [0, weights(8)], 'r');
plot([points(9), points(9)], [0, weights(9)], 'r')
plot([points(10), points(10)], [0, weights(9)], 'r')
plot([points(11), points(11)], [0, weights(10)], 'r')
plot([points(12), points(12)], [0, weights(10)], 'r')
plot([points(13), points(13)], [0, weights(11)], 'r')
plot([points(14), points(14)], [0, weights(11)], 'r')
legend([p1,p2,p3], 'Rosen/Hong 2 point', 'Hong 3 point', 'Zhao 7 point')
hold off
title('X = Rayleigh(1)')
x = sort(x);
y = chi2pdf([0:0.1:9], 2);
figure()
plot(0:0.1:9,y)
title('Y = Chi-Square(2)')


a =5;
b =2;
x = betarnd(a,b,1,n);
yfun = @(x)x+2;
pdf_x = betapdf(x_2, a, b);
[x_params, analytic, sampled, pem_direct, pem_iter, hong, hongs, zhao, points, weights] = pem(yfun, x, x_2, pdf_x, 'B', a, b);
disp('X = Beta(5,2), Y = X + 2')
weights(1:4)
points(1:4)
table(x_params, analytic, sampled, pem_direct, pem_iter, hong, hongs, zhao, 'RowNames', rows, 'VariableNames', cols)
x = sort(x);
figure()
hold on
plot(x, betapdf(x,a,b))
p1 = plot([points(1), points(1)], [0, weights(1)], 'g');
plot([points(2), points(2)], [0, weights(2)], 'g')

p2 = plot([points(5), points(5)], [0, weights(5)], 'b');
plot([points(6), points(6)], [0, weights(6)], 'b')
plot([points(7), points(7)], [0, weights(7)], 'b')

p3 = plot([points(8), points(8)], [0, weights(8)], 'r');
plot([points(9), points(9)], [0, weights(9)], 'r')
plot([points(10), points(10)], [0, weights(9)], 'r')
plot([points(11), points(11)], [0, weights(10)], 'r')
plot([points(12), points(12)], [0, weights(10)], 'r')
plot([points(13), points(13)], [0, weights(11)], 'r')
plot([points(14), points(14)], [0, weights(11)], 'r')
legend([p1,p2,p3], 'Rosen/Hong 2 point', 'Hong 3 point', 'Zhao 7 point')
hold off
title('X = Beta(5,2)')
y = betapdf(x, a, b);
figure()
plot(x + 2,y)
title('Beta displace by +2')

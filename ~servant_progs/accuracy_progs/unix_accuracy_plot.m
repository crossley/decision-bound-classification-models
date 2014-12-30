% Script for making contour plot of accuracy for unix rule

xc_vector = 0:100;
xc_vector = xc_vector+0.01;

[accuracy, max_accuracy, xc_max] = comp_unix_acc_fcn(xc_vector);

figure
plot(xc_vector, accuracy)



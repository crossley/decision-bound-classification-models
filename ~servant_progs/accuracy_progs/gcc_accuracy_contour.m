% Script for making contour plot of accuracy for gcc rule

xc_vector = 0:100;
yc_vector = 0:100;

[accuracy, max_accuracy, xc_max, yc_max] = comp_gcc_acc_fcn(xc_vector, yc_vector);

[X,Y] = meshgrid(xc_vector, yc_vector);

figure
contour(X,Y,accuracy,20)



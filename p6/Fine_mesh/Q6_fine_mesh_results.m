clc 
clear
close all


stresses = readtable('fpprts.dat');
stresses.Properties.VariableNames = ["Output","V1","Element ID","V2","sigma_zz","sigma_xx","sigma_yy","tau_xy"];
x_axis_coordinates = readtable('Element_coordinates_GP_x_axis.xlsx');
y_axis_coordinates = readtable('Element_coordinates_GP_y_axis.xlsx');

%% Theoretical expressions
sigma_0 = 1;
a = 1/2;
sigma_r = @(r,theta) sigma_0/2*(1-a^2/r^2)+sigma_0/2*(1+3*a^4/r^4-4*a^2/r^2)*cos(2*theta);
sigma_theta = @(r,theta) sigma_0/2*(1+a^2/r^2)-sigma_0/2*(1+3*a^4/r^4)*cos(2*theta);


%% Stress variation along x-axis
stresses_x_axis = stresses(1:28,:);

figure()
plot(x_axis_coordinates.x,stresses_x_axis.sigma_yy,'LineWidth',2)
hold on
fplot(@(r) sigma_theta(r,0),[0.5 5],'--r','LineWidth',2)
xlabel('x-axis [mm]')
ylabel('$$\sigma_{yy}$$ [MPa]','Interpreter','latex')
legend('Fine Mesh FEA','Analytic Result')
grid on

%% Stress variation along y-axis
stresses_y_axis = stresses(29:end,:);

figure()
plot(y_axis_coordinates.y,stresses_y_axis.sigma_xx,'LineWidth',2)
hold on
fplot(@(r) sigma_theta(r,pi/2),[0.5 2.5],'--r','LineWidth',2)
xlabel('y-axis [mm]')
ylabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
legend('Fine Mesh FEA','Analytic Result')
grid on

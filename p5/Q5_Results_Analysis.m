clear
close all
clc

% Data
P = 2;
L = 8;
c = 1;
I = 2*c^3/3;
E = 1000;
nu = 0.25;

% Theoretical Results
y_vec = linspace(-c,c,100);

%% x = 0.25
% sigma_xx
x = 0.25;
sigma_xx_1_theo = -P.*y_vec.*(L-x)/I;

% tau_xy
tau_xy_1_theo = P*(c^2-y_vec.^2)/(2*I);
%% x = 3.75
x = 3.75;
sigma_xx_2_theo = -P.*y_vec.*(L-x)/I;

% tau_xy
tau_xy_2_theo = P*(c^2-y_vec.^2)/(2*I);

% vertical tip deflection
v_tip_theo = P/(6*E/(1-nu^2)*I)*((L-L)^3-L^3+L*((4+5*nu/(1-nu))*c^2 + 3*L^2));



%% Part a)
%% x = 0.25
% 2-elem mesh
sigma_xx_1_2_elem = [11.0023491391793 -11.0023491391793];
y_vec_2_elem = [-0.5 0.5];
tau_xx_1_2_elem = [1 1];

% 4-elem mesh
sigma_xx_1_4_elem = [16.8784918971520 5.83445985852038 -5.83445985852000 -16.8784918971517];
y_vec_4_elem = [-0.75 -0.25 0.25 0.75];
tau_xx_1_4_elem = [0.501647925313714  1.49835207468567 1.49835207468587 0.501647925313594];

% 8-elem mesh
sigma_xx_1_8_elem = [19.6139822860202 14.4158555323487 8.68676463671505 ...
    2.89975810915660 -2.89975810915546 -8.68676463671481 ... 
-14.4158555323484 -19.6139822860192];
y_vec_8_elem = [-0.875 -0.625 -0.375 -0.125 0.125 0.375 0.625 0.875];
tau_xx_1_8_elem = [0.146681834722543 0.929003556995897 1.36057270582963 ...
    1.56294190244888 1.56294190244943 1.36057270582981 0.929003556995886 ...
    0.146681834722190];

% 2-elem quad mesh
%% x = 0.5
sigma_xx_1_2_elem_quad = [11.2505390715217 -11.2505390715219];
y_vec_2_elem = [-0.5 0.5];
tau_xx_1_2_elem_quad = [0.971862509986697 0.971862509986681];

figure()
plot(sigma_xx_1_2_elem_quad,y_vec_2_elem,'LineWidth',2)
hold on
plot(-P.*y_vec.*(L-0.5)/I,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem quad mesh','Analytical Result')
%title('$$\sigma_{xx}$$ for x = 0.50','Interpreter','latex')

figure()
plot(tau_xx_1_2_elem_quad,y_vec_2_elem,'LineWidth',2)
hold on
plot(tau_xy_1_theo,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem quad mesh','Analytical Result')
%title('$$\tau_{xy}$$ for x = 0.50','Interpreter','latex')

%% x = 3.5
sigma_xx_2_2_elem_quad = [6.75003472294435 -6.75003472294448];
y_vec_2_elem = [-0.5 0.5];
tau_xx_2_2_elem_quad = [0.969493080411543 0.969493080411554];

figure()
plot(sigma_xx_2_2_elem_quad,y_vec_2_elem,'LineWidth',2)
hold on
plot(-P.*y_vec.*(L-3.5)/I,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem quad mesh','Analytical Result')
%title('$$\sigma_{xx}$$ for x = 0.50','Interpreter','latex')

figure()
plot(tau_xx_1_2_elem_quad,y_vec_2_elem,'LineWidth',2)
hold on
plot(tau_xy_2_theo,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem quad mesh','Analytical Result')
%title('$$\tau_{xy}$$ for x = 0.50','Interpreter','latex')

% Plots

figure()
plot(sigma_xx_1_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(sigma_xx_1_4_elem,y_vec_4_elem,'LineWidth',2)
plot(sigma_xx_1_8_elem,y_vec_8_elem,'LineWidth',2)
plot(sigma_xx_1_theo,y_vec,'LineWidth',2)
ylabel('y [mm]')
xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','Analytical Result')
%title('$$\sigma_{xx}$$ for x = 0.25','Interpreter','latex')

figure()
plot(tau_xx_1_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(tau_xx_1_4_elem,y_vec_4_elem,'LineWidth',2)
plot(tau_xx_1_8_elem,y_vec_8_elem,'LineWidth',2)
plot(tau_xy_1_theo,y_vec,'LineWidth',2)
ylabel('y [mm]')
xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','Analytical Result')
%title('$$\tau_{xy}$$ for x = 0.25','Interpreter','latex')









%% x = 3.75
% 2-elem mesh
sigma_xx_2_2_elem = [6.04444444444843 -6.04444444444836];
y_vec_2_elem = [-0.5 0.5];
tau_xx_2_2_elem = [1 1];

% 4-elem mesh
sigma_xx_2_4_elem = [9.27272727365251 3.09090928036503 -3.09090928036473 -9.27272727365225];
y_vec_4_elem = [-0.75 -0.25 0.25 0.75];
tau_xx_2_4_elem = [0.636363313506699  1.36363668649246 1.36363668649245 0.636363313506688];

% 8-elem mesh
sigma_xx_2_8_elem = [10.8778222907826 7.76987486349608 4.66192747563307 ...
    1.55397637196408 -1.55397637196355 -4.66192747563254 ... 
-7.76987486349554 -10.8778222907820];
y_vec_8_elem = [-0.875 -0.625 -0.375 -0.125 0.125 0.375 0.625 0.875];
tau_xx_2_8_elem = [0.359926300847468 0.908387749499634 1.27403083293760 ...
    1.45685511671253 1.45685511671256 1.27403083293770 0.908387749499701 ...
    0.359926300847563];


figure()
plot(sigma_xx_2_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(sigma_xx_2_4_elem,y_vec_4_elem,'LineWidth',2)
plot(sigma_xx_2_8_elem,y_vec_8_elem,'LineWidth',2)
plot(sigma_xx_2_theo,y_vec,'LineWidth',2)
ylabel('y [mm]')
xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','Analytical Result')
%title('$$\sigma_{xx}$$ for x = 3.75','Interpreter','latex')

figure()
plot(tau_xx_2_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(tau_xx_2_4_elem,y_vec_4_elem,'LineWidth',2)
plot(tau_xx_2_8_elem,y_vec_8_elem,'LineWidth',2)
plot(tau_xy_2_theo,y_vec,'LineWidth',2)
ylabel('y [mm]')
xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','Analytical Result')
%title('$$\tau_{xy}$$ for x = 3.75','Interpreter','latex')



%% vertical tip deflection
v_tip_2_elem = 0.477167;
v_tip_4_elem = 0.484520;
v_tip_8_elem = 0.485953;
v_tip_2_elem_quad = 5.01235E-01;


%% Part b)
% Different BCs
%% x = 0.25
% 8-elem mesh
sigma_xx_1_8_elem_bcs = [20.3973365121318 13.6821467658986 8.01365018518896 ...
    2.64434774379137 -2.64434774379062 -8.01365018518822 ... 
-13.6821467658979 -20.3973365121310];
y_vec_8_elem = [-0.875 -0.625 -0.375 -0.125 0.125 0.375 0.625 0.875];
tau_xx_1_8_elem_bcs = [1.62270331750376 1.15123062512135 0.717779329656609 ...
    0.507486727713927 0.507486727713914 0.717779329656567 1.15123062512127 ...
    1.62270331750367];

%% x = 3.75
sigma_xx_2_8_elem_bcs = [10.8778086965578 7.76988994502078 4.66193038807201 ...
    1.55399146903345 -1.55399146903279 -4.66193038807134 ... 
-7.76988994502015 -10.8778086965570];
y_vec_8_elem = [-0.875 -0.625 -0.375 -0.125 0.125 0.375 0.625 0.875];
tau_xx_2_8_elem_bcs = [0.359918784365637 0.908386473732820 1.27403037656907 ...
    1.45686436532959 1.45686436532958 1.27403037656910 0.908386473732770 ...
    0.359918784365676];

figure()
plot(sigma_xx_1_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(sigma_xx_1_4_elem,y_vec_4_elem,'LineWidth',2)
plot(sigma_xx_1_8_elem,y_vec_8_elem,'LineWidth',2)
plot(sigma_xx_1_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
plot(sigma_xx_1_theo,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','8-elem mesh w/ diff BCs','Analytical Result')


% plot(sigma_xx_1_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
% hold on 
% plot(sigma_xx_1_theo,y_vec,'LineWidth',2)
% ylabel('y [mm]')
% xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
% grid on
% legend('8-elem mesh w/ diff BCs','Analytical Result')
% title('$$\sigma_{xx}$$ for x = 0.25','Interpreter','latex')

figure()
plot(tau_xx_1_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(tau_xx_1_4_elem,y_vec_4_elem,'LineWidth',2)
plot(tau_xx_1_8_elem,y_vec_8_elem,'LineWidth',2)
plot(tau_xx_1_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
plot(tau_xy_1_theo,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','8-elem mesh w/ diff BCs','Analytical Result')
% plot(tau_xx_1_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
% hold on 
% plot(tau_xy_1_theo,y_vec,'LineWidth',2)
% ylabel('y [mm]')
% xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
% grid on
% legend('8-elem mesh w/ diff BCs','Analytical Result')
% title('$$\tau_{xy}$$ for x = 0.25','Interpreter','latex')

figure()
plot(sigma_xx_2_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(sigma_xx_2_4_elem,y_vec_4_elem,'LineWidth',2)
plot(sigma_xx_2_8_elem,y_vec_8_elem,'LineWidth',2)
plot(sigma_xx_2_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
plot(sigma_xx_2_theo,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','8-elem mesh w/ diff BCs','Analytical Result')

% figure()
% plot(sigma_xx_2_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
% hold on 
% plot(sigma_xx_2_theo,y_vec,'LineWidth',2)
% ylabel('y [mm]')
% xlabel('$$\sigma_{xx}$$ [MPa]','Interpreter','latex')
% grid on
% legend('8-elem mesh w/ diff BCs','Analytical Result')
% title('$$\sigma_{xx}$$ for x = 0.50','Interpreter','latex')

figure()
plot(tau_xx_2_2_elem,y_vec_2_elem,'LineWidth',2)
hold on 
plot(tau_xx_2_4_elem,y_vec_4_elem,'LineWidth',2)
plot(tau_xx_2_8_elem,y_vec_8_elem,'LineWidth',2)
plot(tau_xx_2_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
plot(tau_xy_2_theo,y_vec,'--r','LineWidth',2)
ylabel('y [mm]')
xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
grid on
legend('2-elem mesh','4-elem mesh','8-elem mesh','8-elem mesh w/ diff BCs','Analytical Result')

% figure()
% plot(tau_xx_2_8_elem_bcs,y_vec_8_elem,'LineWidth',2)
% hold on 
% plot(tau_xy_2_theo,y_vec,'LineWidth',2)
% ylabel('y [mm]')
% xlabel('$$\tau_{xy}$$ [MPa]','Interpreter','latex')
% grid on
% legend('8-elem mesh w/ diff BCs','Analytical Result')
% title('$$\tau_{xy}$$ for x = 0.50','Interpreter','latex')


%% vertical deflection
v_tip_8_elem_bcs = 0.486724;

%% Part c)
% vertical deflection
v_tip_distorted = 0.48125; 

% Output of the deflection
fprintf('The theoretical tip deflection is %.4f mm.\n',v_tip_theo)
fprintf('The 2-elem mesh gives a tip deflection of %.4f mm with an error of %.4f %%.\n',v_tip_2_elem, 100*abs((v_tip_2_elem-v_tip_theo)/v_tip_theo));
fprintf('The 4-elem mesh gives a tip deflection of %.4f mm with an error of %.4f %%.\n',v_tip_4_elem, 100*abs((v_tip_4_elem-v_tip_theo)/v_tip_theo));
fprintf('The 8-elem mesh gives a tip deflection of %.4f mm with an error of %.4f %%.\n',v_tip_8_elem, 100*abs((v_tip_8_elem-v_tip_theo)/v_tip_theo));
fprintf('The 2-elem quad mesh gives a tip deflection of %.4f mm with an error of %.4f %%.\n',v_tip_2_elem_quad, 100*abs((v_tip_2_elem_quad-v_tip_theo)/v_tip_theo));
fprintf('The 8-elem mesh with different BCs gives a tip deflection of %.4f mm with an error of %.4f %%.\n',v_tip_8_elem_bcs, 100*abs((v_tip_8_elem_bcs-v_tip_theo)/v_tip_theo));
fprintf('The mesh with distorted element gives a tip deflection of %.4f mm with an error of %.4f %%.\n',v_tip_distorted, 100*abs((v_tip_distorted-v_tip_theo)/v_tip_theo));


%% BONUS
nu = 0.499;
% vertical tip deflection
v_tip_theo = P/(6*E/(1-nu^2)*I)*((L-L)^3-L^3+L*((4+5*nu/(1-nu))*c^2 + 3*L^2))


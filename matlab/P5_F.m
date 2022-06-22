close all
clear all
clc

% GRIPEN-2

% Aircraft Data
h=1000 ; M=0.40; aa0=2.46; gg0=15; u0=261.4; flaps=0;
th0=48; de0=-7.44; da0=0.58; dr0=-0.01;
Teng=0.50; demax=[-22 28]; damax=[-18 18]; drmax=[-23 23]; flapmax=[-40 40];
m=9211; Ix=1783099; Iy=65106; Iz=1734269; Ixz=1763;
S=25.55; b=8.382; c=2.235; aamax=19.69;

% Derivatives - Longitudinal
xu=-0.0384; xw=0.0455; zu=-0.1521; zw=-1.1253; zwp=-0.0045; zq=-1.6353;
mu=0.0000; mw= -0.0803; mq=-1.0798; mwp=-0.0333;
xde=0.000; zde=-0.012; mde=-0.020; xdsp=-1.677; zdsp=0.000; mdsp=0.000;
xdt=8.760; zdt=0.000; mdt=0.000;

% Derivatives - Lateral
ybb=-0.0689; lbb=-0.0696; nbb=0.3144; yp=0.0007; lp=-0.8419; np=0.0019;
yr=0.0018; lr=0.0084; nr=-0.0256;
lda=-0.702; nda=-0.052; ydr=-0.005; ldr=-0.002; ndr=-0.013;

% Auxiliar
g=9.81;
deg=pi/180;
knot2ms=0.514444444;

% Pre calculations
u0=u0*knot2ms;
w0=u0*aa0*deg;
tt0=aa0+gg0;
ctt0=cos(tt0*deg);
stt0=sin(tt0*deg);
ttt0=tan(tt0*deg);
l_bb=lbb+Ixz/Ix*nbb;
l_p=lp+Ixz/Ix*np;
l_r=lr+Ixz/Ix*nr;
l_da=lda+Ixz/Ix*nda;
l_dr=ldr+Ixz/Ix*ndr;
n_bb=nbb+Ixz/Iz*lbb;
n_p=np+Ixz/Iz*lp;
n_r=nr+Ixz/Iz*lr;
n_da=nda+Ixz/Iz*lda;
n_dr=ndr+Ixz/Iz*ldr;

% Lateral Dynamic Equation
A_AA=[ybb  yp+w0/u0 yr-1 g*ctt0/u0;
      l_bb l_p   l_r   0;
      n_bb n_p   n_r   0;
      0    1     ttt0  0];
%damp(A_AA);

B=[0    ydr;
   l_da l_dr;
   n_da n_dr;
   0    0];

% Exit Equation
C = eye(4);
D = zeros(4,2);

% Bryson Method - Final
bb_B   = 0.18;
p_B    = 1;  % 0.5 - experimentar esta com atuadores          
r_B    = 0.4;              
phi_B  = 1;
dda_B  = 2;  % 1 - experimentar esta com atuadores      
ddr_B  = 14;

% Cost Matrix with Bryson Method
Q = diag([1/(bb_B)^2  1/(p_B)^2 1/(r_B)^2 1/(phi_B)^2]);
R = diag([1/(dda_B)^2 1/(ddr_B)^2]);

% LQR
K_lqr = lqr(A_AA,B,Q,R);  % gain matrix
% Closed Loop Dynamic Matrix (6 state Matrix)
A_f = A_AA - B*K_lqr;
damp(A_f)

% Compute F matrix
CC = [1 0 0 0; 0 0 0 1]; % beta and phi
dcgain = dcgain(A_f,B,C,D); % closed loop gain
F = (CC*dcgain)^(-1); % F matrix

% Compute K_psi gain
% G_phi_phiref=tf(ss(A_AA-B*K_lqr,B*F*[0;1],[0 0 0 1],0)); % G = phi/phiref
% s=tf('s');
% G_psi_phiref = G_phi_phiref/(s*(u0/g)); % G = psi/phiref
% rlocus(G_psi_phiref)
% grid
K_psi = 2.49;

% LOS - characteristics
global actualpath Rmin distmin points;
actualpath = 1;
Rmin = 3500; % m
distmin = 1400; % m 800
points=         [0         1.0000
                 1.0000    3.0000
                 1.0000    5.5000
                 1.0000    8.0000
                 3.0000    8.0000
                 3.0000    5.5000
                 3.0000    3.0000
                 5.0000    3.0000
                 5.0000    5.5000
                 5.0000    8.0000
                 7.0000    8.0000
                 7.0000    5.5000
                 7.0000    3.0000
                % 4.2500    0.6000
                 1.5000   -1.8000
                 0.8000   -1.8000
                 0.0000   -1.0000
                 0         1.0000] .* Rmin; % m


% Run simulink
T_final = 1050; %s 1075
bb_stept = 10; %s
bb_ref = 0; %deg
out = sim('P5_F_sim',T_final);

% Plots beta & beta_ref
figure()
plot(out.beta.time,out.beta.data,'r','Linewidth',1.2)
hold on
plot(out.beta_ref.time,out.beta_ref.data,'b','Linewidth',1.2)
legend('\beta','\beta_{ref}','Location','NorthEast')
grid on
xlabel('time [s]')
ylabel('Deg')

% Plots phi & phi_ref
figure()
plot(out.phi_ref.time,out.phi_ref.data,'b','Linewidth',1.2)
hold on
plot(out.phi.time,out.phi.data,'r','Linewidth',1.2)
legend('\phi_{ref}','\phi','Location','NorthEast')
grid on
xlabel('time [s]')
ylabel('Deg')

% Plots delta_a & delta_r
figure()
plot(out.da.time,out.da.data,'b','Linewidth',1.2)
hold on
plot(out.dr.time,out.dr.data,'r','Linewidth',1.2)
legend('\delta_a','\delta_r','Location','NorthEast')
grid on
xlabel('time [s]')
ylabel('Deg')

% Plots p & r
figure()
plot(out.p.time,out.p.data,'b','Linewidth',1.2)
hold on
plot(out.r.time,out.r.data,'r','Linewidth',1.2)
legend('p','r','Location','NorthEast')
grid on
xlabel('time [s]')
ylabel('Deg/s')

% Plots psi & psi_ref
figure()
plot(out.psi.time,out.psi.data,'r','Linewidth',1.2)
hold on
plot(out.psi_ref.time,out.psi_ref.data,'--b','Linewidth',1.2)
legend('\psi','\psi_{ref}','Location','NorthEast')
grid on
xlabel('time [s]')
ylabel('Deg')

% Plots path & path_ref
path = path_ref();
figure()
plot(points(:,1),points(:,2),'*b')
hold on
plot(out.pos.data(:,1),out.pos.data(:,2),'-r','Linewidth',1.2)
plot(path(1,:),path(2,:),'--b','Linewidth',1.2)
xlabel('East [m]')
ylabel('North [m]')
legend('path_{ref} points','path','path_{ref}','Location','SouthEast')


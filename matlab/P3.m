close all
clear all
%clc

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

% Lateral Dinamic Equation
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

% Integrative States for beta and phi
A2_AA = [A_AA zeros(4,2);
    1 0 0 0 0 0;
    0 0 0 1 0 0];
B2 = [B;
    0 0;
    0 0];

% Bryson Method
bb_B   = 0.5 * deg;
p_B    = 0.4 * deg;              
r_B    = 0.18 * deg;              
phi_B  = 0.7 * deg;
ibb_B  = 0.3 * bb_B;
iphi_B = 0.3  * phi_B;
% dimunuir r_B aumenta amortecimento e dimunui fn de RH quase nao afeta R,E
% aumentar ddr_B aumenta amortecimento e fn de RH quase nao afeta R,E

dda_B  = 10 * deg;     
ddr_B  = 20 * deg;

% Q = diag([1/(bb_B)^2  1/(p_B)^2 1/(r_B)^2 1/(phi_B)^2]);
% R = diag([1/(dda_B)^2 1/(ddr_B)^2]);
Q = diag([1/(bb_B)^2  1/(p_B)^2 1/(r_B)^2 1/(phi_B)^2 1/(ibb_B)^2 1/(iphi_B)^2]);
R = diag([1/(dda_B)^2 1/(ddr_B)^2]);

% K_lqr = lqr(A_AA,B,Q,R); 
% A_f = A_AA - B*K_lqr;
% damp(A_f)
K_lqr = lqr(A2_AA,B2,Q,R); 
A2_f = A2_AA - B2*K_lqr;
damp(A2_f)

K_pr = K_lqr(:,2:3);
K_bbphi = [K_lqr(:,1) K_lqr(:,4)];
K_int = K_lqr(:,5:6);

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

B=[0    ydr;
   l_da l_dr;
   n_da n_dr;
   0    0];

damp(A_AA);

% Exit Equation
C = eye(4);
D = zeros(4,2);

% LGR
% LGR - r/dr
figure()
rlocus(A_AA,B(:,2),-[0 0 1 0],0)
sgrid
hold on
kr_dr = 115;  % gain
% plot new poles
p_R1 = rlocus(A_AA,B(:,2),-[0 0 1 0],0,kr_dr);
plot(p_R1,'+k','LineWidth',1.2)
% new dynamic matrix
A_R1 = A_AA + B(:,2)*[0 0 kr_dr 0];
damp(A_R1)

% LGR - bb/dr
figure()
rlocus(A_R1,B(:,2),[1 0 0 0],-0)
sgrid
hold on
kbb_dr = 39;  % gain
% plot new poles
p_R2 = rlocus(A_R1,B(:,2),[1 0 0 0],-0,kbb_dr);
plot(p_R2,'+k','LineWidth',1.2)
% new dynamic matrix
A_R2 = A_R1 - B(:,2)*[kbb_dr 0 0 0];
damp(A_R2)

% LGR - p/da
figure()
rlocus(A_R2,B(:,1),-[0 1 0 0],-0)
sgrid
hold on
kp_da = 1;  % gain 
% plot new poles
p_R3 = rlocus(A_R2,B(:,1),-[0 1 0 0],-0,kp_da);
plot(p_R3,'+k','LineWidth',1.2)
% new dynamic matrix
A_R3 = A_R2 + B(:,1)*[0 kp_da 0 0];
damp(A_R3)


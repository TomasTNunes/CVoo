%-- BAE146 : flight condition : 1

%-- SAE SISO

%Limpa tudo
clear;

pi = 3.14159;
deg = pi/180;
g = 9.81; % m/s^2 
h = 0; %m
M = 0.25; % Mach
aa0 = 5.27 * deg; %em rad 
gg0 = 10 * deg; % em rad
tt0 = (aa0 + gg0); % em rad (cálculo feito por nós)
u0 = 84.99; % em m/s -- //165.2 em kn 
flaps = 12*deg; % em rad 

%throttle : 
th0 = 113; %(%)

de0 = 0.00; % rad
da0 = 0.00; % rad
dr0 = 0.00; % rad

Teng =5.00; % s
demax = +28 * deg; % em rad
demin = -21 * deg; % em rad
damax = 17 * deg; % em rad
drmax = 23 * deg; % em rad
flapmax = 40 * deg; % em rad

%inertial data :
m = 43091; % kg
Ix =16576051; % kg .m ^2
Iy =1926217; % kg .m ^2
Iz =2427183; % kg .m ^2
Ixz =1763; % kg .m ^2

%wing data : 
S = 85.84; % m ^2
b = 26.213; % m
c = 3.277; % m
aamax = 14.10 * deg; % rad

%derivatives ( no units or SI units ):

% estado longitudinal
xu = -0.0310;
xw = 0.0784;
zu = -0.2288;
zw = -0.6659;
zwp = 0.0205;
zq = -3.3631;
mu = 0.0000;
mw = -0.0382;
mq = -1.1856;
mwp = -0.0005;

% estado lateral
ybb = -0.0538; % ybb = yv/u0
yv = ybb; % para colocar na matriz
lbb = -0.0508; % lbb = lv/u0
lv = lbb; % para colocar na matriz
nbb = 0.7144; % nbb = nv/u0
nv = nbb; % para colocar na matriz
yp = 0.0011;
lp = -0.0615;
np = -0.0743;
yr = 0.0069;
lr = -0.0046;
nr = -0.1366;

l_lv = lv + (Ixz/Ix)*nv; % L'v
n_lv = nv + (Ixz/Iz)*lv; % N'v
l_lp = lp + (Ixz/Ix)*np; % L'p
n_lp = np + (Ixz/Iz)*lp; % N'p
l_lr = lr + (Ixz/Ix)*nr; % L'r
n_lr = nr + (Ixz/Iz)*lr; % N'r

% entradas longitudinais
xde = 0.000;
zde = -5.385;
mde = -1.727;
xdf = -0.854;
zdf = -10.597;
mdf = 0.008;
xdt = 3.147;
zdt = 0.000;
mdt = 0.000;

% entradas laterais 
Yda = 0; % assumido
Lda = -0.155;
Nda = 0.138;
Ydr = -0.030;
Ldr = -0.012;
Ndr = -0.489;


% Matrizes do Espçao de Estados

A = [yv (yp+aa0) (yr - 1) (g*cos(tt0)/u0); l_lv l_lp l_lr 0; n_lv n_lp n_lr 0; 0 1 tan(tt0) 0]; 

B = [Yda Ydr ; Lda Ldr ; Nda Ndr ; 0 0];

C = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];

D = [0 0; 0 0; 0 0; 0 0];



%Realimentação do r para o rudder
rlocus(A, B(:,2), -[0 0 1 0], 0);
sgrid;
kr = 1.83;

%Realimentação do p para o aileron
rlocus(A+B(:,2)*[0 0 kr 0], B(:,1), -[0 1 0 0], 0);
sgrid;
kp = 4.7;

damp(A+B(:,2)*[0 0 kr 0]+B(:,1)*[0 kp 0 0]);

open('B_SISO_Simulink');
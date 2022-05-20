% ___________________________________________________
%| Secção 5:                                         |
%|        Aterragem ILS ALTERNATIVA                  |
%|___________________________________________________|
% 
% Apresenta-se agora a segunda alternativa de aterragem, em que os 
% requisitos de estabilidade são cumpridos, mas a modelação dos sensores
% não é incluída.
%
% Nota: antes de executar o presente ficheiro, realizar a instrução 'clear',
% seguida da execucao de 'dados_sensores.m'
%
% A simulacao associada, 'AterragemILS2.m, é executada apos este ficheiro
% ser executado corretamente

%__________________________________________________________________
% Saltar para linha 83
%
% conversoes de unidades
kn=0.514444;   % knot to m/s  
deg=pi/180;    % degree º to rad

% constantes (unidades SI)
h=100;
M=0.12; 
U0=79.2*kn;
aa0=7.17*deg;
gg0=0;
tt0=gg0+aa0;
th0 = 0.73;
Ixx=3094421;
Iyy=4301930;
Izz=4524396;
Ixz=1356;

Ip=Ixz/Ixx;
Ir=Ixz/Izz;
g=9.814;
W0=U0*aa0;
bb0=0; r0=0; p0=0; phi0=0;
da0=0; dr0=0;

v_wind = 10;
absangle_wind = 20; % absolute angle; angle_wind = -20º

% derivadas (unidades SI)
Ybb=-0.0276;
Yp=0.0008;
Yr=0.0001;
YdA=0;          % nao especificado, assumido como 0
YdR=-0.008;
Lbb=-0.0486;
Lp=-0.3703;
Lr=-0.0689;
LdA=-0.194;
LdR=-0.010;
Nbb=0.1060;
Np=0.0118;
Nr=-0.0751;
NdA=0.005;
NdR=-0.058;

%Definiçao dos modelos dinamicos

Lbb2=Lbb+Ip*Nbb;
Lp2=Lp+Ip*Np;
Lr2=Lr+Ip*Nr;
LdA2=LdA+Ip*NdA;
LdR2=LdR+Ip*NdR;

Nbb2=Nbb+Ir*Lbb;
Np2=Np+Ir*Lp;
Nr2=Nr+Ir*Lr;
NdA2=NdA+Ir*LdA;
NdR2=NdR+Ir*LdR;

% Superficies de controlo
dda_max = 20 * deg; %[rad]
ddr_max = 20 * deg; %[rad]
%____________________________________________________________________

% Formulação em espaço de estados 

A=[Ybb  Yp+W0/U0      Yr-1  (g*cos(tt0))/U0  
   Lbb2    Lp2          Lr2             0            
   Nbb2    Np2          Nr2             0  
   0       1           tan(tt0)         0 ];

B=[YdA      YdR      
   LdA2     LdR2
   NdA2     NdR2
   0        0];

C=eye(4);

D=zeros(4,2);

 % Inclusao de um quinto estado na dinamica: ângulo de rumo
 
A3 =[Ybb      Yp+W0/U0     Yr-1                g*cos(tt0)/U0            0              
    Lbb2     Lp2            Lr2                         0               0           
    Nbb2     Np2            Nr2                         0               0           
    0        1              tan(tt0)                    0               0           
    Ybb     Yp+W0/U0      (Yr-1)+1/cos(tt0)   g*cos(tt0)/U0             0];

B3=[YdA      YdR        
   LdA2     LdR2
   NdA2     NdR2
   0        0
   YdA     YdR];

C3 = eye(5,5);
D3 = zeros(5,2);

% Estabilizacao com LQR

bb_B   = 5.1* deg;
p_B    = 1* deg;
r_B    = 2* deg; 
phi_B  = 13.5 * deg;
lambda_B  = 1.1 * deg; 

dda_B  = 30 * deg; 
ddr_B  = 32* deg;

Q = diag([1/(bb_B)^2  1/(p_B)^2 1/(r_B)^2 1/(phi_B)^2 1/(lambda_B)^2]);
R = diag([1/(dda_B)^2 1/(ddr_B)^2]);

k3 = lqr(A3,B3,Q,R);                    
anel_fechado3 = A3-B3*k3;
%damp(A3-B3*k3) 
% Desta vez, a aeronave foi estabilizada com sucesso

k_interior = [k3(:,1:4)];
k_lambda = k3(2,5);

G=-([1 0 0 0 0 ; 0 0 0 1 0 ])*inv(anel_fechado3)*B3;
F=pinv(G);

                            

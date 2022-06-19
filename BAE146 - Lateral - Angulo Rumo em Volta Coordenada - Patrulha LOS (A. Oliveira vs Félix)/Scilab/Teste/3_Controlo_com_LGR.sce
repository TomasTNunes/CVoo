//-- BAE146 : flight condition : 1

//Limpa tudo
clear;


pi = 3.14159;
deg = pi/180;
g = 9.81; // m/s^2 
h = 0; //m
M = 0.25; // Mach
aa0 = 5.27 * deg; //em rad 
gg0 = 10 * deg; // em rad
tt0 = (aa0 + gg0); // em rad (cálculo feito por nós)
u0 = 84.99; // em m/s -- //165.2 em kn 
flaps = 12*deg; // em rad 

//throttle : 
th0 = 113; //(%) 

de0 = 0.00; // rad
da0 = 0.00; // rad
dr0 = 0.00; // rad

Teng =5.00; // s
demax = +28 * deg; // em rad ???
demin = -21 * deg; // em rad ???
damax = 17 * deg; // em rad
drmax = 23 * deg; // em rad
flapmax = 40 * deg; // em rad

//inertial data :
m = 43091; // kg
Ix =16576051; // kg .m ^2
Iy =1926217; // kg .m ^2
Iz =2427183; // kg .m ^2
Ixz =1763; // kg .m ^2

//wing data : 
S = 85.84; // m ^2
b = 26.213; // m
c = 3.277; // m
aamax = 14.10 * deg; // rad

//derivatives ( no units or SI units ):

// estado longitudinal
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

// estado lateral
ybb = -0.0538; // ybb = yv/u0
yv = ybb; //para coloar na matriz
lbb = -0.0508; // lbb = lv/u0
lv = lbb; //para coloar na matriz
nbb = 0.7144; // nbb = nv/u0
nv = nbb; //para coloar na matriz
yp = 0.0011;
lp = -0.0615;
np = -0.0743;
yr = 0.0069;
lr = -0.0046;
nr = -0.1366;

l_lv = lv + (Ixz/Ix)*nv // L'v
n_lv = nv + (Ixz/Iz)*lv //N'v
l_lp = lp + (Ixz/Ix)*np //L'p
n_lp = np + (Ixz/Iz)*lp //N'p
l_lr = lr + (Ixz/Ix)*nr //L'r
n_lr = nr + (Ixz/Iz)*lr //N'r

// entradas longitudinais
xde = 0.000;
zde = -5.385;
mde = -1.727;
xdf = -0.854;
zdf = -10.597;
mdf = 0.008;
xdt = 3.147;
zdt = 0.000;
mdt = 0.000;

//entradas laterais 
Yda = 0; //assumido
Lda = -0.155;
Nda = 0.138;
Ydr = -0.030;
Ldr = -0.012;
Ndr = -0.489;


function rt(Ft) //Função para chamar rootlocus com grelha
    evans(Ft);
    sgrid("red");
endfunction


A = [yv (yp+aa0) (yr - 1) (g*cos(tt0)/u0) 0; l_lv l_lp l_lr 0 0; n_lv n_lp n_lr 0 0; 0 1 tan(tt0) 0 0; 0 0 1/cos(tt0) 0 0]; 

B = [Yda Ydr ; Lda Ldr ; Nda Ndr ; 0 0; 0 0];

D = [0 0; 0 0; 0 0; 0 0; 0 0]

x0 = [0; 0; 0; 0; 0]

C = diag([1,1,1,1,1])

s = %s

//psi = xx;

Kp = 4.72;
Kr = 1.825;
Kpsi = 1;
Kphi = 1;
psi_ref = 0;

K = [0 Kp 0 -Kphi Kphi*Kpsi*-1); 0 0 Kr 0 0];

esp_est = syslin('c', A + B*K,B,C,D,x0); 

pol = spec(esp_est(2))
[wn, zz] = damp(pol)

// Voltar a determinar os níveis de qualidade de voo === [Classe III] e [Categoria A]

//modo espiral
T2 = log(2) / wn(5,1) // tem que ser maior que 12s
re = real(pol(5,1))
if re < 0
    disp("Espiral estável")
    disp(pol(5,1))
elseif T2 < 12 then
    disp("T2 Espiral menor que 12");
    disp(T2);
end

//modo rolamento
Tp = 1 / wn(3,1) // tem que ser menor 1.4
if 1.4 < Tp then
    disp("Tp Rolamento maior que 1.4");
    disp(Tp);
end

//modo rolamento holandes
zz_rh = zz(2,1) // maior do que 0.19 (imposto pelo projeto maior que 0.6)
if zz_rh < 0.6 then
    disp("Amortecimento RH menor que 0.6");
    disp(zz_rh);
end

wn_rh = wn(2,1) // maior do que 0.5
if wn_rh < 0.5 then
    disp("Frequencia RH menor que 0.5");
    disp(wn_rh);
end

zz_wn_rh = zz_rh * wn_rh // maior do que 0.35
if zz_wn_rh < 0.35 then
    disp("Freq*Amort RH menor que 0.35");
    disp(zz_wn_rh);
end

//Mostrar todos os gráficos de todas as saídas em resposta de todas as entradas

be=B(:,1); // 1 - aileron 2 - rudder
c=C(1,:); // bb p r phi
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(0); 
t1=0:0.001:100;
gs1=csim('impulse',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de bb", "fontsize", 2, "color","blue");
xtitle('Resposta LQR de bb a um step no aileron');
/*
c=C(2,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(1); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de p", "fontsize", 2, "color","blue");
xtitle('Resposta LQR da razão de rolamento a um step no aileron');

c=C(3,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(2); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de r", "fontsize", 2, "color","blue");
xtitle('Resposta LQR razão de guinada a um step no aileron');

c=C(4,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(3); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de phi", "fontsize", 2, "color","blue");
xtitle('Resposta LQR do ângulo phi a um step no aileron');

c=C(5,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(4); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de phi", "fontsize", 2, "color","blue");
xtitle('Resposta LQR do ângulo phi a um step no rudder');

be=B(:,2); // 1 - aileron 2 - rudder
c=C(1,:); // bb p r phi
sysaaq=syslin('c',A + B*K,be,c,0,x0);

figure(5); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de bb", "fontsize", 2, "color","blue");
xtitle('Resposta LQR de bb a um step no rudder');

c=C(2,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(6); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de p", "fontsize", 2, "color","blue");
xtitle('Resposta LQR da razão de rolamento a um step no rudder');

c=C(3,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(7); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de r", "fontsize", 2, "color","blue");
xtitle('Resposta LQR razão de guinada a um step no rudder');

c=C(4,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(8); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de phi", "fontsize", 2, "color","blue");
xtitle('Resposta LQR do ângulo phi a um step no rudder');

c=C(5,:);
sysaaq=syslin('c',A + B*K,be,c,0,x0);
figure(9); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step' == rudder ou aileron
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de phi", "fontsize", 2, "color","blue");
xtitle('Resposta LQR do ângulo phi a um step no rudder');

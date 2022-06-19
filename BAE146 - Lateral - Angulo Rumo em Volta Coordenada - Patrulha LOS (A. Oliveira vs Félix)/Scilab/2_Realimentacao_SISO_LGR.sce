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
yv = ybb; //para colocar na matriz
lbb = -0.0508; // lbb = lv/u0
lv = lbb; //para colocar na matriz
nbb = 0.7144; // nbb = nv/u0
nv = nbb; //para colocar na matriz
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


A = [yv (yp+aa0) (yr - 1) (g*cos(tt0)/u0); l_lv l_lp l_lr 0; n_lv n_lp n_lr 0; 0 1 tan(tt0) 0]; 

B = [Yda Ydr ; Lda Ldr ; Nda Ndr ; 0 0];

D = [0 0; 0 0; 0 0; 0 0]

x0 = [0; 0; 0; 0]

C = diag([1,1,1,1]);


function rt(Ft) //Função para chamar rootlocus com grelha
    evans(Ft);
    sgrid("red");
endfunction


//Modelo completo 
// == LGR da realimentação da razão de guinada para a defleção do rudder no Rolamento Holandês
// == LGR da realimentação da razão de rolamento para a defleção do aileron no Rolamento

//  Tp = 1/ wn = 1.4 == wn = 0.71 == polo rolamento = -0.71 == Ganho visto no RL p/delta_a ---- rt(G(2,1))
//  zz_rh tem que ser maior que 0.6 == ganho visto no RL r/delta_r ---- rt(G(3,2))
//olhando para o root locus:
//       ~se quisermos aumentar o zz de rh, é impossível com realimentação neg.
//       ~se quisermos aumentar o wn de r, é impossível com realimentação neg.
// ===== SÃO NECESSÁRIAS DUAS REALIMENTAÇÕES POSITIVAS

Kr = 1.825 //visto no root locus rt(G(3,2)) em anel aberto

Kp = 4.72 // visto no root locus rt(G(2,1)) em anel aberto

// realimentação da razão de guinada ([0 0 1 0]) para a defleção do rudder (B(:,2)) com ganho Kr:
// +
// realimentação da razão de rolamento ([0 1 0 0]) para a defleção do a (B(:,2)) com ganho Kp
esp_est = syslin('c', A + B(:,2)*Kr*C(3,:) + B(:,1)*Kp*C(2,:),B,C,D,x0); 

G = -ss2tf(esp_est); //menos(-) para fazer ver feedback positivo no root locus
pol = spec(esp_est(2))
[wn, zz] = damp(pol)

// Voltar a determinar os níveis de qualidade de voo === [Classe III] e [Categoria A]

//modo espiral
T2 = log(2) / wn(4,1) // tem que ser maior que 12s
re = real(pol(4,1))
if re < 0
    disp("Espiral estável")
    disp(pol(4,1))
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

/*
Com estas realimentações:
+Todos os modos estão no nível 1 das qualidaades de voo.
*/

/*
be=B(:,2); // 1 - aileron 2 - rudder
c=[0 0 1 0]; // bb p r phi
sysaaq=syslin('c',A + B(:,2)*Kr*C(3,:) + B(:,1)*Kp*C(2,:),be,c,0,x0);
//resposta a um degrau
figure(0); 
t1=0:0.001:100;
gs1=csim('step',t1,sysaaq); // 'impuls', 'step'
plot2d(t1,gs1);
xlabel("tempo(seg)", "fontsize", 2, "color","blue"); 
ylabel("Amplitude de r", "fontsize", 2, "color","blue");
xtitle('Resposta do Anel Aberto da razão de picada a um step');
*/

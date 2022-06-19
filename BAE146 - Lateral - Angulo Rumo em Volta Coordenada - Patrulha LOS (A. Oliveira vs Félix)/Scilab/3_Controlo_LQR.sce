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


A = [yv (yp+aa0) (yr - 1) (g*cos(tt0)/u0) 0; l_lv l_lp l_lr 0 0; n_lv n_lp n_lr 0 0; 0 1 tan(tt0) 0 0; 0 0 1/cos(tt0) 0 0]; 

B = [Yda Ydr ; Lda Ldr ; Nda Ndr ; 0 0; 0 0];

D = [0 0; 0 0; 0 0; 0 0; 0 0]

x0 = [0; 0; 0; 0; 0]

C = diag([1,1,1,1,1]);

//Método de Bryson (valores máximos arbitrados por nós)

max_bb = 2*deg; // bb = asin(vento_longitunial/u0): maximizado == vento_longitunial/u0 maximo == asin(10/84.99)  
max_p = 1*deg; // 0.025
max_r = 5*deg; // 0.01
max_phi = 5*deg; //  dado (max volta coordeanada)
max_psi = 25*deg; //prof indicou 5x-10x o phi
max_da = 25*deg; // prof indicou 10% do maximo absoluto fornecido colocamos 60% (volta coordenada usa ailerons) 0.35
max_dr = 25*deg; // prof indicou 10% retirado do maximo absoluto fornecido 0.17

Q = diag([1/(max_bb^2);1/(max_p^2);1/(max_r^2);1/(max_phi^2);1/(max_psi^2)]) //max esperado para: bb; p; r; phi
R = diag([1/(max_da^2);1/(max_dr^2)]) //max esperado para: delta_a; delta_r

esp_est = syslin('c',A,B,C,D,x0); 

//Método LQR para determinar ganhos
K = lqr(esp_est,Q,R);

esp_est = syslin('c', A + B*K,B,C,D,x0); 

pol = spec(esp_est(2))
[wn, zz] = damp(pol)
disp(pol)

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
Tp = 1 / wn(1,1) // tem que ser menor 1.4
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

//Variavel auxiliar
global ponto;
ponto = 1;
//Funcao para calcular o angulo necessario para ir ter ao ponto pretendido
function angulopsi = atan2(y,x) // E,N
    if (0 < x && 0 < y) then
        angulopsi = atan(y/x);
    elseif (0 < x && y < 0) then
        angulopsi = atan(y/x) + 2*%pi;
    elseif (x < 0 && 0 <= y) then
        angulopsi = atan(y/x) + %pi;
    elseif (x < 0 && y < 0) then
        angulopsi = atan(y/x) + %pi;
    elseif (x == 0 && y < 0) then
        angulopsi = 3*%pi/2;
    elseif (x == 0 && 0 < y) then
        angulopsi = %pi/2;
    elseif (y == 0 && x < 0) then
        angulopsi = %pi;
    elseif (y == 0 && 0 < x) then
        angulopsi = 0;
    elseif (x == 0 && y == 0) then
        angulopsi = 0;
    end
endfunction

// Funcao para definir a trajetoria
function angulo_psi = path(posicao_x,posicao_y)
    global ponto;
    //Definicao dos pontos do trajeto
    x_path = [0;1181;8900;24463];
    y_path = [0;3529;18245;18245];
    //Diferenca entre os posicao atual e ponto do trajeto
    delta_x = x_path(ponto,1) - posicao_x;
    delta_y = y_path(ponto,1) - posicao_y;
    //Calculo da distancia entre a posicao atual e o ponto do trajeto
    distance = sqrt(delta_x^2 + delta_y^2);
    //Avaliacao da posicao atual relativamente ao ponto do trajeto
    if  distance < 300 then //quando se aproximar do ponto, muda para o proximo ponto do trajeto
        ponto = ponto+1;
        if ponto == 4 then
            ponto = 1;
        end
        delta_x = x_path(ponto,1) - posicao_x;
        delta_y = y_path(ponto,1) - posicao_y;
        distance = sqrt(delta_x^2 + delta_y^2);
    end
    angulo_psi = atan2(delta_x,delta_y)
endfunction

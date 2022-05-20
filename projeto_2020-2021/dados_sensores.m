% Conversoes de unidades
kn=0.514444;        % knot to m/s  
deg=pi/180;         % degree to rad

% Frequencia de amostragem
f = 40; %[Hz]
    
% Atuadores [SI]
act_t = 0.1; %[s]
act_vmax = 1; %[rad/s]

% S1: Angulos Aerodinamicos [SI]
s1_inmax=25*deg; %[rad]
s1_inmin=-25*deg; %[rad]
s1_outmin=0; %[V]
s1_outmax=5; %[V]
s1_k = (s1_outmax-s1_outmin)/(s1_inmax-s1_inmin);
s1_offset = s1_outmax-s1_k*s1_inmax;
s1_t_const=0.01; %[s]
s1_noise = 0.005; %[V]
s1_power=s1_noise^2/f;

% S2: Razoes angulares [SI]
s2_inmax=50*deg; %[rad/s]
s2_inmin=-50*deg; %[rad/s]
s2_outmin=-3; %[V]
s2_outmax=3; %[V]
s2_k = (s2_outmax-s2_outmin)/(s2_inmax-s2_inmin);
s2_offset = s2_outmax-s2_k*s2_inmax;
s2_noise = 0.002; %[V]
s2_power=s2_noise^2/f;

% S3: Angulo de rolamento [SI]
s3_inmax=90*deg; %[rad]
s3_inmin=-90*deg; %[rad]
s3_outmin=0; %[V]
s3_outmax=28; %[V]
s3_k = (s3_outmax-s3_outmin)/(s3_inmax-s3_inmin);
s3_offset = s3_outmax-s3_k*s3_inmax;
s3_noise = 0.25*deg; %[rad]
s3_power = s3_noise^2/f;

% S4: Angulo de guinada [SI]
s4_inmax=360*deg; %[rad]
s4_inmin=0*deg; %[rad]
s4_outmin=0; %[V]
s4_outmax=28; %[V]
s4_k = (s4_outmax-s4_outmin)/(s4_inmax-s4_inmin);
s4_offset = s4_outmax-s4_k*s4_inmax;
s4_noise = 1.5*deg; %[rad]
s4_power=s4_noise^2/f;

% S5: ILS [SI]              
s5_sensibilidade= 3.63*10^-6; %[A/º]
s5_max = 150*10^-6; %[A]



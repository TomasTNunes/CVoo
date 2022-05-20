% Atencao: Correr esta ficheiro apenas depois de correr EKC39.m
%
% ___________________________________________________
%| Secção 2:                                         |
%|        Aumento de estabilidade                    |
%|___________________________________________________|
%
% Neste ficheiro acompanha-se o processo de aumento de estabilidade da
% aeronave

%ROOT_LOCUS
%yaw_damper
sys_sae_a= ss( A, B(:,2), -[0 0 1 0], 0);
%rlocus(sys_sae_a), grid, axis equal
k_a=[0 0 0 0;
     0 0 4 0];%4

A_a_iteration= A +B*k_a;
%damp(A_a_iteration)


sys_sae_b = ss(A_a_iteration, B(:,1), -[0 1 0 0],0);
%rlocus(sys_sae_b), grid, axis equal
k_b=[0 2 0 0;%2
     0 0 0 0];
A_b_iteration= A_a_iteration + B*k_b;
%damp(A_b_iteration)


sys_sae_c= ss( A_b_iteration, B(:, 1), +[0 0 1 0], 0);
%rlocus(sys_sae_c), grid, axis equal
k_c=[0 0 1 0;
     0 0 0 0];
A_c_iteration=A_b_iteration + B*k_c;
%damp(A_c_iteration)
 
k=k_a+k_b;


%LQR USADO PARA MOSTRAR QUE A ESTABILIZACAO DA AERONAVE E POSSIVEL

%Q=diag([  (2*deg)^-2    (1*deg)^-2    (1*deg)^-2    (10*deg)^-2   ]);
%R=diag([ (2*deg)^-2     (3*deg)^-2]);

Q=diag([  (2.5*deg)^-2    (1*deg)^-2    (1*deg)^-2    (10*deg)^-2   ]);
R=diag([ (10*deg)^-2     (10*deg)^-2]);
k_lqr_1=lqr(A,B,Q,R);
damp(A-B*k_lqr_1);
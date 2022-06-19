function out=E_rumo_ref(posicao)

    global rumoactual; %E           N
    pontos=1.0e+003 *[  0           0.850
                        1.700 		2.550    
                        1.700		6.800
                        3.400       6.800
                        3.400		2.550 
                        5.100       2.550
                        5.100       6.800
                        6.800       6.800
                        6.800       2.550
                        2.000      -1.500
                       -0.200      -1.500
                        0           0.900];
    
   
    rumo=[pontos(rumoactual,1) pontos(rumoactual,2)];
    dist=sqrt((rumo(1)-posicao(1))^2+(rumo(2)-posicao(2))^2);

    if(dist<=1400)
        rumoactual=rumoactual+1;
        if(rumoactual==12)
        end
        rumo=[pontos(rumoactual,1) pontos(rumoactual,2)];
    end
        angulo=atan2((rumo(1)-posicao(1)), (rumo(2)-posicao(2)));
     
    %para corre��o da descontinuidade do atan2
    if(angulo<=0)
        out=angulo+2*pi;
    end
    
    if(angulo>0)
        out=angulo;
    end 

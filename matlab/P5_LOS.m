function out=P5_LOS(posicao)
    global rumoactual Rmin distmin pontos;
    
    rumo = [pontos(rumoactual,1) pontos(rumoactual,2)];
    dist = sqrt((rumo(1)-posicao(1))^2+(rumo(2)-posicao(2))^2);

    if dist <= distmin
        rumoactual = rumoactual + 1;
        if rumoactual == 12
        end
        rumo=[pontos(rumoactual,1) pontos(rumoactual,2)];
    end
    angulo = atan2((rumo(1) - posicao(1)), (rumo(2) - posicao(2)));
     
    %para correção da descontinuidade do atan2
    if(angulo<=0)
        out = angulo + 2*pi;
    end
    
    if(angulo>0)
        out = angulo;
    end 

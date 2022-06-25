function out=P5_LOS(pos)
% Computes psi ref
    global actualpath Rmin distmin points;
    
    rumo = [points(actualpath,1) points(actualpath,2)]; % LOS point
    dist = sqrt((rumo(1)-pos(1))^2+(rumo(2)-pos(2))^2);

    if dist <= distmin
        actualpath = actualpath + 1
        if actualpath == length(points)
            actualpath = 1;
        end
        rumo=[points(actualpath,1) points(actualpath,2)]; % LOS point
    end

    angulo = atan2((rumo(1) - pos(1)), (rumo(2) - pos(2))); % psi ref
     
    % Correction of atan2 between 0º and 360º
    if(angulo<=0)
        out = angulo + 2*pi;
    end
    
    if(angulo>0)
        out = angulo;
    end 

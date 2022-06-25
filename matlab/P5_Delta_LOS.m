function out=P5_LOS(pos)
    global actualpath Rmin distmin points ref_path Delta LOS;
    
    rumo = [points(actualpath+1,1) points(actualpath+1,2)];
    dist = sqrt((rumo(1)-pos(1))^2+(rumo(2)-pos(2))^2);

    if dist <= distmin
        actualpath = actualpath + 1
        if actualpath == length(points)
            actualpath = 1;
        end
    end
    
    POS = [pos(1) pos(2)];
    P = @(k) points(actualpath,:) + ref_path(actualpath,:) .* k;
    dotp = @(k) dot( P(k)-points(actualpath,:) , POS-P(k) );
    options = optimset('Display','off');
    K = fsolve(dotp,1,options);
    PLOS = P(K) + ref_path(actualpath,:) / norm(ref_path(actualpath,:)) * Delta;
    if norm(PLOS-points(actualpath,:)) > norm(ref_path(actualpath,:))
        PLOS = points(actualpath+1,:);
    end
    angulo = atan2((PLOS(1) - pos(1)), (PLOS(2) - pos(2)));
    LOS=[LOS; PLOS];
     
    %para correção da descontinuidade do atan2
    if(angulo<=0)
        out = angulo + 2*pi;
    end
    
    if(angulo>0)
        out = angulo;
    end 

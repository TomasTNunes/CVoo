function path = path_ref_Delta()
% Computes Reference path vectors
global points Rmin
pointsE = points(:,1);
pointsN = points(:,2);
pathE = [pointsE(1) pointsE(2) pointsE(3)];
pathN = [pointsN(1) pointsN(2) pointsN(3)];

c1 = @(x) pointsN(3) + sqrt(Rmin^2 - (x - (pointsE(3) + pointsE(4))/2).^2);
pathE = [pathE linspace(pointsE(3),pointsE(4),100)];
pathN = [pathN c1(linspace(pointsE(3),pointsE(4),100))];

pathE = [pathE pointsE(5)];
pathN = [pathN pointsN(5)];

c2 = @(x) pointsN(5) - sqrt(Rmin^2 - (x - (pointsE(5) + pointsE(6))/2).^2);
pathE = [pathE linspace(pointsE(5),pointsE(6),100)];
pathN = [pathN c2(linspace(pointsE(5),pointsE(6),100))];

pathE = [pathE pointsE(7)];
pathN = [pathN pointsN(7)];

c3 = @(x) pointsN(7) + sqrt(Rmin^2 - (x - (pointsE(7) + pointsE(8))/2).^2);
pathE = [pathE linspace(pointsE(7),pointsE(8),100)];
pathN = [pathN c3(linspace(pointsE(7),pointsE(8),100))];

pathE = [pathE pointsE(9) pointsE(10) pointsE(11) pointsE(12)];
pathN = [pathN pointsN(9) pointsN(10) pointsN(11) 0];

path = [pathE; pathN];
end
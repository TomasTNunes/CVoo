function path = path_ref()

global points Rmin
pointsE = points(:,1);
pointsN = points(:,2);
pathE = [pointsE(1) pointsE(2) pointsE(4)];
pathN = [pointsN(1) pointsN(2) pointsN(4)];

c1 = @(x) pointsN(4) + sqrt(Rmin^2 - (x - (pointsE(4) + pointsE(5))/2).^2);
pathE = [pathE linspace(pointsE(4),pointsE(5),100)];
pathN = [pathN c1(linspace(pointsE(4),pointsE(5),100))];

pathE = [pathE pointsE(7)];
pathN = [pathN pointsN(7)];

c2 = @(x) pointsN(7) - sqrt(Rmin^2 - (x - (pointsE(7) + pointsE(8))/2).^2);
pathE = [pathE linspace(pointsE(7),pointsE(8),100)];
pathN = [pathN c2(linspace(pointsE(7),pointsE(8),100))];

pathE = [pathE pointsE(10)];
pathN = [pathN pointsN(10)];

c3 = @(x) pointsN(4) + sqrt(Rmin^2 - (x - (pointsE(10) + pointsE(11))/2).^2);
pathE = [pathE linspace(pointsE(10),pointsE(11),100)];
pathN = [pathN c3(linspace(pointsE(10),pointsE(11),100))];

pathE = [pathE pointsE(13) pointsE(14) pointsE(15) pointsE(16) pointsE(17)];
pathN = [pathN pointsN(13) pointsN(14) pointsN(15) pointsN(16) 471];

path = [pathE; pathN];
end
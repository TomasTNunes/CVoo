function out = correct_err(err)
% Angle between -180º and 180º
if err > pi
    out = err - 2*pi;
elseif err < -pi
    out = err + 2*pi;
else
    out = err;
end

end
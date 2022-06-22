function out = correct_err(err)

if err > pi
    out = err - 2*pi;
elseif err < -pi
    out = err + 2*pi;
else
    out = err;
end

end
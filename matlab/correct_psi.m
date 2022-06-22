function out = correct_psi(psi)

if psi > 2*pi
    out = psi - 2*pi;
elseif psi < 0
    out = psi + 2*pi;
else
    out = psi;
end

end
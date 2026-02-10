function f = fermi_dirac(E, Ef, kT)
%FERMI_DIRAC Fermi-Dirac occupation for energies in eV.
%
% kT <= 0 时退化为 0K 阶跃分布。

    if kT <= 0
        f = double(E < Ef);
    else
        f = 1 ./ (1 + exp((E - Ef)/kT));
    end
end


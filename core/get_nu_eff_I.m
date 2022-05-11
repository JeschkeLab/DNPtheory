function [nu_eff_I,kI] = get_nu_eff_I(nu_I,nu_m)

kI = round(nu_I./nu_m);
nu_eff_I = nu_I - kI.*nu_m;

end


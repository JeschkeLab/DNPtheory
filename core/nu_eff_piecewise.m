function [nu_eff,z_eff,nu_m] = nu_eff_piecewise(tp_vec,phi_vec,nu1_vec,offset_vec)
%NU_EFF_PIECEWISE [nu_eff] = nueff_piecewise(tp_vec,phi_vec,nu1_vec,offset_vec)
%   OUTPUT: nu_eff: magnitude of the effective field of the sequence
%           z_eff:  direction of the effectice field vec (normalized)
%   INPUT:  tp_vec:  vector of pulses (and delays)
%            phi_vec: phase of pulses in radians
%           nu1_vec: amplitude of pulses (or delays if 0) in linear freqs
%           offset_vec: offset in MHz (can either be a vector or a single number)
         
%check input (could be expanded...)
if length(offset_vec)==1 && length(tp_vec)>1
   offset_vec = offset_vec*ones(size(tp_vec)); 
end

%generate quaternions for each individual pulse/delay
q_pulse = zeros(4,numel(tp_vec));
for it = 1:numel(tp_vec)
    q_pulse(:,it)=quat_rf(nu1_vec(it),phi_vec(it),offset_vec(it),tp_vec(it),'frame') ;
end

%multiply them step by step (note that it is a frame trafo, hence the order)
% multiply them and
% express them as an array of rotation matrices
qtot = [1 0 0 0]';
for it=1:numel(tp_vec)
    qtot = quatmult(qtot,q_pulse(:,it));
end

%get modulation frequency from the length of the sequence
Tmod = sum(tp_vec);
nu_m = 1/Tmod;

%get effective field magnitude and direction
[z_eff,xi_eff] = quat2eff(qtot);
weff = xi_eff/Tmod;
nu_eff = weff/(2*pi);


end


function [quat] = quat_rf(nu1,phi,offset,tp,direction)
% [quat] = quat_rf(w1,phi,Omega,tp) generates a unit quaternion from the
% rf-field defined by the rf amplitude w1 (in agular units), the phase 
% theta (theta=0: along x-axis), and the offset angular frequency Omega.
% Applied for a duration tp. 
%
% direction:
% 'frame' -> frame transformation, inverts the effective field sign
% 'prop' -> propagation
% quat is a unit quaternion, the first entry is cos(angle/2)

w1=2*pi*nu1;
Omega=offset*2*pi;

weff = sqrt(Omega^2+w1^2);
if strcmp(direction,'frame')
   weff = -weff; 
end

beff = weff*tp;

qr = cos(beff/2);
s2 = sin(beff/2);

if Omega~=0
    theta = atan2(w1,Omega);
    
    qi = s2*cos(phi)*sin(theta);
    qj = s2*sin(phi)*sin(theta);
    qk = s2*cos(theta);
else
%     theta = pi/2; %not needed really
    qi= s2*cos(phi);
    qj = s2*sin(phi);
    qk = 0;
end

quat = [qr qi qj qk]';

end

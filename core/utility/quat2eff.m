function [axis,winkel] = quat2eff(quat)
%[axis,winkel] = quat2eff(quat) extracts the axis of rotation axis and the
%the effective rotation angle from a quaternion. the last three entries of
%quat specify the rotation axis.


max_beta_eff = pi;


if size(quat,1)==1
    quat = quat';
end

thresh= 1e-6;

winkel = 2*acos(quat(1));

if abs(winkel)>thresh && (2*pi-abs(winkel))>thresh
    ax = quat(2:4);
    axis = ax/sqrt(ax'*ax);
else
    axis=[0 0 1]' ;
    winkel=0;
end


if winkel>max_beta_eff
    winkel = 2*pi-winkel;
    axis = - axis;
end


if axis(3)==-1
    axis(3)=1;
    winkel=-winkel;
end

end


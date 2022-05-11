function [varargout]=DNP_numerical(tp_vec,phi_vec,nu1_vec,offset_vec,nu_I,rho0_vec,r,Ntheta,Nrounds)

if length(offset_vec)==1 && length(tp_vec)>1
   offset_vec = offset_vec*ones(size(tp_vec)); 
end

Sx = sop([1/2 1/2],'xe');
Sy = sop([1/2 1/2],'ye');
Sz = sop([1/2 1/2],'ze');


Ix = sop([1/2 1/2],'ex');
Iy = sop([1/2 1/2],'ey');
Iz = sop([1/2 1/2],'ez');

natural_constants
T = mu0/(4*pi)*gfree*bmagn*g1H*nmagn*1/(r*1e-10)^3/planck/1e6;


theta_vec = linspace(0,pi/2,Ntheta);
weights=sin(theta_vec); weights=weights/sum(weights);

sig = zeros(1,Nrounds);

for itheta = 1:Ntheta
    
    theta=theta_vec(itheta);
    A=T*(3*cos(theta)^2-1);
    B=3*sin(theta)*cos(theta)*T;
    
    H0 = nu_I*Iz+A*Sz*Iz+B*Sz*Ix;
    
    Uround = eye(size(Sx));
    for ip = 1:numel(tp_vec)
        H1 = offset_vec(ip)*Sz + nu1_vec(ip)*(cos(phi_vec(ip))*Sx + sin(phi_vec(ip))*Sy);
        Uround = expm(-1i*2*pi*(H0+H1)*tp_vec(ip))*Uround;
    end
 
    rho = rho0_vec(1)*Sx+rho0_vec(2)*Sy+rho0_vec(3)*Sz;
    detect = Iz;
    
    for istep = 1:Nrounds
        sig(istep)=sig(istep) + weights(itheta)*real(trace(detect*rho)/trace(detect*detect));
        rho = Uround*rho*Uround';
    end
    
end

sig_num=sig;
t = (0:(Nrounds-1))*(sum(tp_vec));


if nargout==1
    varargout{1}=sig_num;
else 
    varargout{1}=t;
    varargout{2}=sig_num;
end


end


function [ rf, time ] = build_rf(tp,nu1,phi,dt )
%BUILD_RF
%rf complex represenation of nu1(t)



%check whether input makes sense
if length(tp)==length(nu1) && length(tp)==length(phi)
    N=length(tp);
else
    warning('You have to specify tp, nu1 and phi for each pulse!')
end


%make pulses multiple of basic timestep for the pulses
n_pulse=zeros(1,N);
for n=1:N
    n_pulse(n) = round(tp(n)/dt);
    tp_new = n_pulse(n)*dt;
    if(abs(tp_new-tp(n))>10e-10)
        warning(['set tp(' num2str(n) ') to ' num2str(tp_new) '.']);
        tp(n)=tp_new;
    end
end

%
% build rf
%

rf=zeros(1,sum(n_pulse)); %% yes this plus is needed (because we start at t=0)

start=1;
for n=1:N %loop over all pulses
    for ii=start:(start+n_pulse(n)-1)
        rf(ii)=nu1(n)*exp(1i*phi(n));
    end
    start=ii+1;
end



time=0:dt:(length(rf)-1)*dt;


%     keyboard
end


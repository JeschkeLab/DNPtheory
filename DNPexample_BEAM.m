clear, close all

% add helper functions
addpath(genpath('./core/'))


%% define sequence parameters
nu_I=-14.83;            % Nuclear Zeeman frequency
nu_1=32;                  % Electron Rabi freq
tp1=20*1e-3;            % first pulse length in BEAM
tp2=(1+tp1*(nu_1-abs(nu_I)))/(abs(nu_I)+nu_1);             % second pulse length in BEAM


nu1_vec = nu_1*[1 1];
tp_vec = [tp1 tp2];
phi_vec = [0 180]/180*pi; % phases of the pulses
dt =0.1e-4;             % time step of numerical IFT

offset = 5;

rho0_vec = [-1 0 0]';

%% build rf
[ rf, time ]=build_rf(tp_vec,nu1_vec,phi_vec,dt);

%plotting
h = figure(1);
clf
hold on
plot(time*1e3,real(rf),'b')
plot(time*1e3,imag(rf),'r')
xlabel('t / ns')
ylabel('\nu_1 / MHz')
axis([-1 60 -40 40])
legend('real(rf)','imag(rf)','location','northeast')

% h=make_plot_nice(h);
% exportgraphics(h,'SI_rf.pdf','BackgroundColor','none')

%% IFT, calculate R_control(t)

% generate a vector of offsets the same length as rf
offset=offset*ones(size(rf));

%calculate modulation frequency
wmod=2*pi/(time(end)+dt);
nu_m = wmod/2/pi;


% pre-allocate rotation matrices
R_control = zeros(3,3,numel(time));
R_flipped = zeros(3,3,numel(time));


% build array of rotation quaternions for every step
q_pulse = zeros(4,numel(time));
for it = 1:numel(time)
    q_pulse(:,it)=quat_rf(abs(rf(it)),angle(rf(it)),offset(it),dt,'frame') ;
end


% multiply the quaternionas step-by-step and
% express them as an array of rotation matrices (for plotting only at this stage)
qtot = [1 0 0 0]';
q_control = zeros(4,numel(time));
for it=1:numel(time)
    q_control(:,it)=qtot;
    R_control(:,:,it) = quat2rotmat(qtot);
    qtot = quatmult(qtot,q_pulse(:,it));
end



% Plot IF in rotating frame
IF_plot(time,R_control,2);

% h2=figure(2);
% exportgraphics(h2,'SI_R_control.pdf','BackgroundColor','none')



%% flip each quaternion of q_control (or R_control) 
% such that the effective field is now along the z-axis

%determine the effective field (effective flip angle and direction)
[z_eff,beta_eff] = quat2eff(qtot);

% determine R_flip, or the corresponding quaternion, respectively
% see https://doi.org/10.1063/1.5123046 for details
u_z = [0 0 1]';
q_flip = cross(z_eff,u_z);
qr_flip = 1+u_z'*z_eff;
q_flip = [qr_flip; q_flip ];
q_flip = q_flip/sqrt(sum(q_flip'*q_flip));

% flip all quaternions/rot-matrices
q_flipped= zeros(4,numel(time));
for it=1:numel(time)
    q_flipped(:,it)=quatmult(q_flip,q_control(:,it));
    R_flipped(:,:,it) = quat2rotmat(q_flipped(:,it));
end

%% remove the effective field from the IFT, i.e. calculate R_eff(t)
% determine the effective field frequency
w_eff = beta_eff/(time(end)+dt);
nu_eff = w_eff/(2*pi);

%pre-allocate rotation quaternions/matrices
q_eff = zeros(4,numel(time));
R_eff = zeros(3,3,numel(time));

% calculate R_eff(t)
for it=1:numel(time)
    
    %calculate Rz(-weff*t), 
    beta = -w_eff*time(it);
    q_z = [cos(beta/2) sin(beta/2)*[0 0 1]]';
    
    q_eff(:,it)=quatmult(q_z,q_flipped(:,it));
    R_eff(:,:,it) = quat2rotmat(q_eff(:,it));
end

IF_plot(time,R_eff,3);
% h3=figure(3);
% exportgraphics(h3,'SI_R_eff.pdf','BackgroundColor','none')


%%  calculate Fourier coefficients a^(k)
N=numel(time);
k_vec = ((0:N-1)-fix(N/2));
if mod(N,2)==0
    k_vec = k_vec(2:end);
end

A_k = zeros(3,3,numel(k_vec));


% Fourier transform each element (xx,xy,xz,...zz) of R_eff(t)
for ii=1:3
    for jj=1:3
        y = squeeze(R_eff(ii,jj,:));
        a = fftshift(fft(y))/numel(y);
        if mod(N,2)==0
            A_k(ii,jj,:)=a(2:end);
        else
            A_k(ii,jj,:)=a(1:end);
        end
    end
end


%extract the original z operator (should no oszillate during free evolution)
Az_k = squeeze(A_k(:,3,:));

%% extract the relevant scaling factor

[nu_eff_I,kI] = get_nu_eff_I(nu_I,nu_m);

kmax=max(k_vec);

ax = squeeze(Az_k(1,:));
ay = squeeze(Az_k(2,:));
ap = ax-1i*ay;
am = ax+1i*ay;
az = squeeze(Az_k(3,:));


k0 = kmax+1;
k = kmax+1+kI;
mk = kmax+1-kI;

rel_sign = sign(nu_eff*nu_eff_I);
if rel_sign == 1
    a_eff= -sqrt(ap(mk)*am(k)) ;
else
    a_eff= sqrt(am(mk)*ap(k)) ;
end

proj = z_eff'*rho0_vec;

f_pm = a_eff*proj

%% comparison with numerical simulation

Ntheta = 31; %number of orientations

r=4.5; % e-n distance in Angstrom
Nrounds = 200; %number of repetitions of the DNP element

%caluclate the transfer numerically
[t_num,sig_num]=DNP_numerical(tp_vec,phi_vec,nu1_vec,offset,nu_I,rho0_vec,r,Ntheta,Nrounds);

h4=figure(4);
clf
hold on
plot(t_num,sig_num)

%% analytical calculation
%calculate the mismatch
delta_nu_eff = abs(nu_eff_piecewise(tp_vec,phi_vec,nu1_vec,offset))-abs(nu_eff_I);

%calculate the anisotropy of the hf-coupling
natural_constants
T = mu0/(4*pi)*gfree*bmagn*g1H*nmagn*1/(r*1e-10)^3/planck/1e6;

%generate orientations and corresponding weights
theta_vec = linspace(0,pi/2,Ntheta);
weights=sin(theta_vec); weights=weights/sum(weights);

%preallocate results
t_theo = t_num(1:10:end);
sig_theo = zeros(1,numel(t_theo));

%loop over orientations, sum up results
for itheta = 1:Ntheta
    
    theta=theta_vec(itheta);
    B=3*sin(theta)*cos(theta)*T;
    nu_pm = sqrt(B^2*a_eff^2/4+delta_nu_eff^2); %transfer frequency
    transfer_amp = a_eff^2*B^2/(4*nu_pm^2); %amplitude of transfer
    
    sig_theo = sig_theo+weights(itheta)*sign(a_eff)*proj*transfer_amp*...
                sin(1/2*(2*pi*nu_pm)*t_theo).^2;
end
plot(t_theo,sig_theo,'o')

axis([xlim -0.1 0.8])
xlabel('t / \mus')
ylabel('<I_z>')
legend('num.','theo.','location','northeast')

% h4=make_plot_nice(h4);
% exportgraphics(h4,'theo_vs_num.pdf','BackgroundColor','none')
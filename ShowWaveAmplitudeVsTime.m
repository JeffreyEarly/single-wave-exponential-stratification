file = '/Users/jearly/Documents/ProjectRepositories/single-wave-exponential-stratification/output_181113_decomp.nc';
t = ncread(file, 't');
nT = length(t);

i = 5;
j = 1;
A_plus = double(squeeze( ncread(file,'Ap_realp', [i j 1], [1 1 nT], [1 1 1])  + sqrt(-1)*ncread(file,'Ap_imagp', [i j 1], [1 1 nT], [1 1 1]) ));
A_minus = double(squeeze( ncread(file,'Am_realp', [i j 1], [1 1 nT], [1 1 1])  + sqrt(-1)*ncread(file,'Am_imagp', [i j 1], [1 1 nT], [1 1 1]) ));

A_plus = double(ncread(file,'Ap_realp')) +  + sqrt(-1)*double(ncread(file,'Ap_imagp'));

dx = 31.25;
Lx = 4000;
p = 3;
T_diss = 120;
k_star = power(-T_diss*log(.5)/t(end),1/(2*p))/(2*dx);

nu = (-1^(p+1)) * power(dx/pi,2*p)/T_diss;
lambda = nu*power(2*pi*i/Lx,2*p);
i_star = 64*power(-120/86400*log(0.99),1/(2*p));
m_star = 512*power(-12/86400*log(0.99),1/(2*p));
% abs(A_plus(1))*exp(-lambda*t)

i = 9;
j = 1;
A_plus2 = double(squeeze( ncread(file,'Ap_realp', [i j 1], [1 1 nT], [1 1 1])  + sqrt(-1)*ncread(file,'Ap_imagp', [i j 1], [1 1 nT], [1 1 1]) ));
A_minus2 = double(squeeze( ncread(file,'Am_realp', [i j 1], [1 1 nT], [1 1 1])  + sqrt(-1)*ncread(file,'Am_imagp', [i j 1], [1 1 nT], [1 1 1]) ));

figure
plot(t/3600, abs(A_plus)), ylog, hold on
plot(t/3600, abs(A_plus2))
% title('Amplitude of primary mode', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
% xlabel('time (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
% ylabel('depth integrated energy (m^3/s^2)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

% print('-depsc', sprintf('PrimaryModeAmplitude.eps'))

return;
figure, plot(angle(A_plus))

theta = angle(A_plus);
dtheta = diff(theta);
dtheta(dtheta<0) = dtheta(dtheta<0) + 2*pi;
theta = [theta(1); cumsum(dtheta)];
omega = (theta(100)-theta(1))/(t(100)-t(1));

wavemodel.Omega(i,j,1)/omega
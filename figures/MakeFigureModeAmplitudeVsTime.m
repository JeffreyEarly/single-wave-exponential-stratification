scaleFactor = 1;
LoadFigureDefaults;

file = '/Users/jearly/Documents/ProjectRepositories/single-wave-exponential-stratification/output_181113_decomp.nc';
t = ncread(file, 't');
nT = length(t);

indices = find(abs(A_minus)/max(max(abs(A_plus)))>0.1);

iTime = floor(1/1);
A_plus = double(squeeze( ncread(file,'Ap_realp', [1 1 iTime], [Inf Inf 1], [1 1 1])  + sqrt(-1)*ncread(file,'Ap_imagp', [1 1 iTime], [Inf Inf 1], [1 1 1]) ));
A_minus = double(squeeze( ncread(file,'Am_realp', [1 1 iTime], [Inf Inf 1], [1 1 1])  + sqrt(-1)*ncread(file,'Am_imagp', [1 1 iTime], [Inf Inf 1], [1 1 1]) ));

spectralEnergy = sum(sum(sum(A_plus.*conj(A_plus) + A_minus.*conj(A_minus))));
fprintf('total spectral energy: %f m^3/s\n', spectralEnergy);

dx = 31.25;
p = 3;
T_diss = 120;
k_star = power(-T_diss*log(.5)/t(end),1/(2*p))/(2*dx);
decayIndex = find( wavemodel.k > k_star,1,'first');

dz = 9.7656;
p = 3;
T_diss = 12;
m_star = 512*power(-T_diss*log(.5)/t(end),1/(2*p));


% figure
% subplot(1,2,1)
% pcolor(abs(A_plus)), shading flat
% subplot(1,2,2)
% pcolor(abs(A_minus)), shading flat


normalize = @(x) log(abs(x));

colorlimits = normalize([max(max(abs(A_plus)))/100 max(max(abs(A_plus)))]);

FigureSize = [50 50 figure_width_2col+8 600*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

nHarmonics = 4;
for iHarmonic = 1:nHarmonics
    i = 1 + (iHarmonic-1)*4;
    j = 1;
    Nj = 15;
    A_plus = normalize(double(squeeze( ncread(file,'Ap_realp', [i j 1], [1 Nj nT], [1 1 1])  + sqrt(-1)*ncread(file,'Ap_imagp', [i j 1], [1 Nj nT], [1 1 1]) )));
    A_minus = normalize(double(squeeze( ncread(file,'Am_realp', [i j 1], [1 Nj nT], [1 1 1])  + sqrt(-1)*ncread(file,'Am_imagp', [i j 1], [1 Nj nT], [1 1 1]) )));
    
  
    subplot( nHarmonics,2,2*(iHarmonic-1) + 1)
    pcolor(t/3600,1:Nj,A_plus), shading flat
    caxis(colorlimits)
    if iHarmonic == 1
       title('A_+', 'FontSize', figure_axis_label_size, 'FontName', figure_font); 
    end
    if iHarmonic < nHarmonics
        set(gca, 'XTick', []);
    else
        xlabel('time (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
    end
    if iHarmonic == 1
        ylabel('vertical mode, k=0', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
    elseif iHarmonic == 2
        ylabel('vertical mode, k=k_0', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
    elseif iHarmonic == 3
        ylabel('vertical mode, k=2k_0', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    elseif iHarmonic == 4
        ylabel('vertical mode, k=3k_0', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    end
    
    subplot(nHarmonics,2,2*(iHarmonic-1) + 2)
    pcolor(t/3600,1:Nj,A_minus), shading flat
    set(gca, 'YTick', []);
    caxis(colorlimits)
    if iHarmonic == 1
       title('A_-', 'FontSize', figure_axis_label_size, 'FontName', figure_font); 
    end
    if iHarmonic < nHarmonics
       set(gca, 'XTick', []); 
    else
        xlabel('time (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
    end
end

print('-dpng', '-r300', sprintf('ModeAmplitudesVsTimeNew.png'))

return

i = 5;
j = 1;
A_plus = double(squeeze( ncread(file,'Ap_realp', [i j 1], [1 1 nT], [1 1 1])  + sqrt(-1)*ncread(file,'Ap_imagp', [i j 1], [1 1 nT], [1 1 1]) ));
A_minus = double(squeeze( ncread(file,'Am_realp', [i j 1], [1 1 nT], [1 1 1])  + sqrt(-1)*ncread(file,'Am_imagp', [i j 1], [1 1 nT], [1 1 1]) ));

netcdf.close(ncid);	


figure, plot(real(A_plus)), hold on, plot(abs(A_plus))
figure, plot(angle(A_plus))

theta = angle(A_plus);
dtheta = diff(theta);
dtheta(dtheta<0) = dtheta(dtheta<0) + 2*pi;
theta = [theta(1); cumsum(dtheta)];
omega = (theta(100)-theta(1))/(t(100)-t(1));

wavemodel.Omega(i,j,1)/omega
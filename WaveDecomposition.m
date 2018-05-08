model_output_folder = '/Users/jearly/Documents/ProjectRepositories/single-wave-exponential-stratification/WintersModelRuns/output_180506';

WM = WintersModel(model_output_folder);

[x,y,z,rho_bar,rho_prime] = WM.VariableFieldsFrom2DOutputFileAtIndex(1,'x','y','z','s1_bar','s1');

[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
rhoFunction = @(z) rhoFunc(z) - rhoFunc(max(zIn)) + double(max(rho_bar));

Nx = length(x);
Ny = 1;
Nz = length(z);
dx = (max(x)-min(x))/(Nx-1);
Lx = double(dx*Nx);
Lz = max(zIn)-min(zIn);
latitude = 0;

z = linspace(min(zIn),max(zIn),Nz);
wavemodel = InternalWaveModelArbitraryStratification([Lx, Lx, Lz], [Nx, Ny, Nz], rhoFunction, z, Nz, latitude);

[u,v,rho_prime] = WM.VariableFieldsFrom2DOutputFileAtIndex(1,'u','v','s1');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(0,u,v,rho_prime);
% Cim notes that the vortexHKE plot appears to show peaks along the axis of
% the radius of deformation for each mode.
%
% Check the individual autocorrelation functions leading to the structure
% in the vortex decorrelation time for the nonlinear run.
%
% In the region of wavenumber-mode space where the wave solutions do not
% decorrelated, how does the relative energy differ from a GM spectrum? In
% the decorrelated region, we don't expect this to match.

scaleFactor = 1;
LoadFigureDefaults;

file = '/Users/jearly/Documents/ProjectRepositories/single-wave-exponential-stratification/WintersModelRuns/output_180506';

WM = WintersModel(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Do the wave-vortex decomposition
%
if WM.NumberOf3DOutputFiles > 1
    wavemodel = WM.wavemodel;
    maxFileIndex = WM.NumberOf3DOutputFiles;
    
    % These first two lines actually do the decomposition
    [t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(maxFileIndex,'t','u','v','w','rho_prime');
    wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);
    

elseif WM.NumberOf2DOutputFiles > 1
    [t_max,x,y,z,rho_bar] = WM.VariableFieldsFrom2DOutputFileAtIndex(WM.NumberOf2DOutputFiles,'t','x','y','z','s1_bar');
    
    z = z-(max(z)-min(z));
    [rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
    rhoFunction = @(z) rhoFunc(z) - rhoFunc(max(zIn)) + double(min(rho_bar));
    

    
    Nx = length(x);
    Ny = 1;
    Nz = length(z);
    dx = (max(x)-min(x))/(Nx-1);
    Lx = double(dx*Nx);
    Lz = max(zIn)-min(zIn);
    latitude = 0;
%     wavemodel = InternalWaveModelExponentialStratification([Lx, Lx, Lz], [Nx, Ny, Nz], [5.2e-3 1300], z, floor(Nz/3.9), latitude);
    
    ProfileDeviation = std((rho_bar-wavemodel.internalModes.rhoFunction(z))/max(rho_bar));
    fprintf('Our assumed profile deviates from the recorded profile by 1 part in 10^{%d}\n',round(log10(ProfileDeviation)));

    [t,u,v,w,rho_prime] = WM.VariableFieldsFrom2DOutputFileAtIndex(1,'t','u','v','w','s1');
    wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);

else
    error('No valid files found')
end

waveHKE_plus_full = wavemodel.Ppm_HKE_factor .*  abs(wavemodel.Amp_plus).^2 ;
waveHKE_minus_full = wavemodel.Ppm_HKE_factor .* abs(wavemodel.Amp_minus).^2;
vortexHKE_full = wavemodel.P0_HKE_factor .* abs(wavemodel.B).^2;

totalEnergy = sum(sum(sum(waveHKE_plus_full+waveHKE_minus_full,'omitnan'),'omitnan'),'omitnan')



% waveHKE_full =  ( abs(wavemodel.Amp_plus).^2 + abs(wavemodel.Amp_minus).^2 );
% vortexHKE_full = abs(wavemodel.B).^2;

% Now we need to figure out how to display this information...

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(wavemodel.Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

% Thi is the final output axis for wavenumber
k = reshape(kAxis(1:(length(kAxis)-1)),[],1);

% Mode axis is just what we already have
j = wavemodel.j;

Kh = wavemodel.Kh;
RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);

nK = length(k);
nModes = length(j);

waveHKE_plus = zeros(nK,nModes);
waveHKE_minus = zeros(nK,nModes);
vortexHKE = zeros(nK,nModes);

for iMode = 1:1:nModes
    for iK = 1:1:nK
        indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1)  & ~squeeze(RedundantCoefficients(:,:,1)) );
        for iIndex = 1:length(indicesForK)
            [i,m] = ind2sub([size(Kh,1) size(Kh,2)],indicesForK(iIndex));
            waveHKE_plus(iK,iMode) = waveHKE_plus(iK,iMode) + waveHKE_plus_full(i,m,iMode);
            waveHKE_minus(iK,iMode) = waveHKE_minus(iK,iMode) + waveHKE_minus_full(i,m,iMode);
            vortexHKE(iK,iMode) = vortexHKE(iK,iMode) + vortexHKE_full(i,m,iMode);
        end
    end
end

waveHKE = waveHKE_plus + waveHKE_minus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the damping scales
%

% Kraig Winters uses an e-fold time to set \nu_x in the hypervicous
% operator. We start by computing \nu_x and \nu_z.
dx = wavemodel.x(2)-wavemodel.x(1);
dz = wavemodel.z(2)-wavemodel.z(1);
nu_x = (-1)^(WM.p_x+1)*power(dx/pi,2*WM.p_x) / WM.T_diss_x;
nu_z = (-1)^(WM.p_z+1)*power(dz/pi,2*WM.p_z) / WM.T_diss_z;

nK = length(wavemodel.k)/2 + 1;
k_diss = abs(wavemodel.k(1:nK));
j_diss = 0:max(wavemodel.j); % start at 0, to help with contour drawing
[K,J] = ndgrid(k_diss,j_diss);
M = (2*pi/(length(wavemodel.j)*dz))*J/2;

% scratch
m = (2*pi/(wavemodel.Nz*dz))*(1:wavemodel.Nz)'/2;
filter_z = exp( nu_z*(sqrt(-1)*m).^(2*WM.p_z)*t_max);

k0 = 4;
N = wavemodel.NumberOfWellConditionedModes(k0,1);
varianceLost_F = zeros(N,1);
varianceGained_F = zeros(N,1);
varianceLost_G = zeros(N,1);
varianceGained_G = zeros(N,1);
allIndices = (1:N)';
for j0=1:N
    otherIndices = setdiff(allIndices,j0);
    
    [xbar,f] = CosineTransformForward(wavemodel.z,wavemodel.Sprime(:,j0,k0));
    xj = wavemodel.Sprime(:,:,k0)\CosineTransformBack(f,filter_z.*xbar);
    
    varianceLost_F(j0) =  xj(j0)^2;
    varianceGained_F(otherIndices) = varianceGained_F(otherIndices) + xj(otherIndices).^2;
    
    [xbar,f] = SineTransformForward(wavemodel.z,wavemodel.S(:,j0,k0));
    xj = wavemodel.S(:,:,k0)\SineTransformBack(f,filter_z(2:end-1).*xbar);
    
    varianceLost_G(j0) =  xj(j0)^2;
    varianceGained_G(otherIndices) = varianceGained_G(otherIndices) + xj(otherIndices).^2;
end


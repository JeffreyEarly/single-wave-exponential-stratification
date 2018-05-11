model_output_folder = '/Users/jearly/Documents/ProjectRepositories/single-wave-exponential-stratification/WintersModelRuns/output_180506';
output_directory = '/Users/jearly/Documents/ProjectRepositories/single-wave-exponential-stratification';

[filepath,name,ext] = fileparts(model_output_folder);
outputfile = fullfile(output_directory,strcat(name,'_decomp.nc'));

WM = WintersModel(model_output_folder);

[x,y,z,rho_bar,rho_prime] = WM.VariableFieldsFrom2DOutputFileAtIndex(1,'x','y','z','s1_bar','s1');

z = z-(max(z)-min(z));
[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
rhoFunction = @(z) rhoFunc(z) - rhoFunc(max(zIn)) + double(min(rho_bar));

% figure, plot(rho_bar,z), hold on, plot(rhoFunction(z),z)

Nx = length(x);
Ny = 1;
Nz = length(z);
dx = (max(x)-min(x))/(Nx-1);
Lx = double(dx*Nx);
Lz = max(zIn)-min(zIn);
latitude = 0;

z = linspace(min(zIn),max(zIn),Nz);
if ~exist('wavemodel','var')
    wavemodel = InternalWaveModelArbitraryStratification([Lx, Lx, Lz], [Nx, Ny, Nz], rhoFunction, z, Nz, latitude,'nEVP',2*Nz+16);
end

nFiles = WM.NumberOf2DOutputFiles;
fileIncrements = 1:1:nFiles;

Nk = length(wavemodel.k);
Nj = length(wavemodel.j);
Nt = length(fileIncrements);

precision = 'double';

if strcmp(precision,'single')
    ncPrecision = 'NC_FLOAT';
    setprecision = @(x) single(x);
    bytePerFloat = 4;
else
    ncPrecision = 'NC_DOUBLE';
    setprecision = @(x) double(x);
    bytePerFloat = 8;
end

cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('SHARE'));
cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(outputfile, cmode);

% Define the dimensions
kDimID = netcdf.defDim(ncid, 'k', Nk);
jDimID = netcdf.defDim(ncid, 'j', Nj);
tDimID = netcdf.defDim(ncid, 't', Nt);

% Define the coordinate variables
kVarID = netcdf.defVar(ncid, 'k', ncPrecision, kDimID);
jVarID = netcdf.defVar(ncid, 'j', ncPrecision, jDimID);
tVarID = netcdf.defVar(ncid, 't', ncPrecision, tDimID);
netcdf.putAtt(ncid,kVarID, 'units', 'radians/m');
netcdf.putAtt(ncid,jVarID, 'units', 'mode number');
netcdf.putAtt(ncid,tVarID, 'units', 's');

% Define the wave-vortex variables
ApRealVarID = netcdf.defVar(ncid, 'Ap_realp', ncPrecision, [kDimID,jDimID,tDimID]);
ApImagVarID = netcdf.defVar(ncid, 'Ap_imagp', ncPrecision, [kDimID,jDimID,tDimID]);
AmRealVarID = netcdf.defVar(ncid, 'Am_realp', ncPrecision, [kDimID,jDimID,tDimID]);
AmImagVarID = netcdf.defVar(ncid, 'Am_imagp', ncPrecision, [kDimID,jDimID,tDimID]);

netcdf.putAtt(ncid,ApRealVarID, 'units', 'm/s');
netcdf.putAtt(ncid,ApImagVarID, 'units', 'm/s');
netcdf.putAtt(ncid,AmRealVarID, 'units', 'm/s');
netcdf.putAtt(ncid,AmImagVarID, 'units', 'm/s');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', wavemodel.latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Lz', wavemodel.Lz);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'source-file', model_output_folder);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

% End definition mode
netcdf.endDef(ncid);

% Add the data for the coordinate variables
netcdf.putVar(ncid, kDimID, wavemodel.k);
netcdf.putVar(ncid, jDimID, wavemodel.j);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 6;
totalSize = totalFields*bytePerFloat*length(fileIncrements)*Nk*Nj/1e9;
fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',outputfile,totalSize);

timeScale = 1; %7.1625; %5.2e-3/7.1890160119109623e-4;

startTime = datetime('now');
for iTime = 1:length(fileIncrements)
    if iTime>=2 %|| mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-1);
        timeRemaining = (length(fileIncrements)-iTime+1)*timePerStep;
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(fileIncrements), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'dd:HH:MM:SS')) ;
    end
    if iTime == 2
        % The first time step takes extra long, because we're using a fixed
        % dimension length for time. So, let's reset the clock for
        % subsequent estimatation.
        startTime = datetime('now');
    end
    [t,u,v,rho_prime] = WM.VariableFieldsFrom2DOutputFileAtIndex(fileIncrements(iTime),'t','u','v','s1');
    
    t = t*timeScale;
    
    wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);
    
    netcdf.putVar(ncid, ApRealVarID, [0 0 iTime-1], [Nk Nj 1], real(squeeze(wavemodel.Amp_plus)));
    netcdf.putVar(ncid, ApImagVarID, [0 0 iTime-1], [Nk Nj 1], imag(squeeze(wavemodel.Amp_plus)));
    netcdf.putVar(ncid, AmRealVarID, [0 0 iTime-1], [Nk Nj 1], real(squeeze(wavemodel.Amp_minus)));
    netcdf.putVar(ncid, AmImagVarID, [0 0 iTime-1], [Nk Nj 1], imag(squeeze(wavemodel.Amp_minus)));
    
    netcdf.putVar(ncid, tVarID, iTime-1, 1, t);
end

netcdf.close(ncid);	

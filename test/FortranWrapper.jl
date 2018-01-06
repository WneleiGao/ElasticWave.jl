# number of grid points
nz = 101; nx = 641;

# grid cell size
 deltaZ = 10.0        ;  deltaX = 10.0;
hdeltaZ = deltaZ / 2.0; hdeltaX = deltaX / 2.0;

# thickness of CPML in number of grid points
nPointPml = 10;

# velocity model
vp = 3300.0; vs = vp / sqrt(2.0); rho = 2800.0;

# time step size
dt = 0.001;

# source parameter
f0  = 7.0;
isz = 20; isx = floor(Int64, nx/2)+1;
zsrc = (isz-1)*deltaZ; xsrc = (isx-1)*deltaX;

# receivers
irz = [50, 50]; irx = [100, 200]

# arrays for elastic model parameter
lambda = zeros(nz+2, nx+2);
mu     = zeros(nz+2, nx+2);
rho    = zeros(nz+2, nx+2);

# arrays for wave field variables
vz = zeros(nz＋2, nx＋2);
vx = zeros(nz＋2, nx＋2);
sigmazz = zeros(nz＋2, nx＋2);
sigmaxx = zeros(nz＋2, nx＋2);
sigmazx = zeros(nz＋2, nx＋2);

# half grid elastic model parameters
hlambdaz    = 0.0
hmuz        = 0.0
hmux        = 0.0
hlambda2muz = 0.0
hrhozx      = 0.0

# parameter for CPML profile
npower   = 2.0
kmax     = 1.0
alphaMax = pi * f0

# damping profile, need to be improved  !!!
 kz=ones(nz);  dz=zeros(nz);  alphaz=zeros(nz);  az=zeros(nz);  bz=zeros(nz);
hkz=ones(nz); hdz=zeros(nz); halphaz=zeros(nz); haz=zeros(nz); hbz=zeros(nz);

 kx=zeros(nx);  dx=zeros(nx);  alphax=zeros(nx);  ax=zeros(nx);  bx=zeros(nx);
hkx=zeros(nx); hdx=zeros(nx); halphax=zeros(nx); hax=zeros(nx); hbx=zeros(nx);

# specify the obsorbing boundary, origin specify the computation domain, z
reflectCoef=1e-3;
d0Z = -(nPower + 1) * vp * log(reflectCoef) / (2.0 * thickPmlZ);
d0X = -(nPower + 1) * vp * log(reflectCoef) / (2.0 * thickPmlX);
locationZ=0.0; locationX=0.0; abscissa=0.0; normalizedAbscissa=0.0;
thickPmlZ = nPointPml * deltaZ; originTop = thickPmlZ; originBottom=(nz-1)*deltaZ - thickPmlZ;

# vertical damping profile
for i = 1 : nz
    locationZ = deltaZ * (i-1)

    # upper profile at grid point
    abscissa  = originTop - locationZ
    if abscissa > 0.0
       normalizedAbscissa = abscissa / thickPmlZ
       dz[i] = d0Z * normalizedAbscissa^nPower
       kz[i] = 1.0
       alphaz[i] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end

    # upper profile at half grid point
    abscissa = originTop - (locationZ + hdeltaZ)
    if abscissa > 0.0
       normalizedAbscissa = abscissa / thickPmlZ
       hdz[i] = d0Z * normalizedAbscissa^npower
       hkz[i] = 1.0
       halphaz[i] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end

    # bottom profile at grid point
    abscissa = locationZ - originBottom
    if abscissa > 0.0
       normalizedAbscissa = abscissa / thickPmlZ
       dz[i] = d0Z * normalizedAbscissa^nPower
       kz[i] = 1.0
       alphaz[i] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end
    # bottom profile at half grid point
    abscissa = (locationZ + hdeltaZ) - originBottom
    if abscissa > 0.0
       normalizedAbscissa = abscissa / thickPmlZ
       hdz[i] = d0Z * normalizedAbscissa^npower
       hkz[i] = 1.0
       halphaz[i] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end

    bz[i]


end




thickPmlX = nPointPml * deltaX; originLeft= thickPmlX; originRight =(nx-1)*deltaX - thickPmlX;





# memory variables for CPML, it can be improved to save memory !!!
mdvzdz = zeros(nz+2, nx+2)
mdvzdx = zeros(nz+2, nx+2)
mdvxdz = zeros(nz+2, nx+2)
mdvxdx = zeros(nz+2, nx+2)
mdsigmazzdz = zeros(nz+2, nx+2)
mdsigmaxxdx = zeros(nz+2, nx+2)
mdsigmazxdz = zeros(nz+2, nx+2)
mdsigmazxdx = zeros(nz+2, nx+2)

# temperal variables to save spatial derivative
vdvzdz = 0.0
vdvzdx = 0.0
vdvxdz = 0.0
vdvxdx = 0.0
vdsigmazzdz = 0.0
vdsigmaxxdx = 0.0
vdsigmazxdz = 0.0
vdsigmazxdx = 0.0


# ==============================================================================
#                              initialization
# ==============================================================================














# crazy

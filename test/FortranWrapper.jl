# number of grid points
nz = 101; nx = 641;

# grid cell size
 deltaZ = 10.0        ;  deltaX = 10.0;
hdeltaZ = deltaZ / 2.0; hdeltaX = deltaX / 2.0;

# time step size
dt = 0.001; nstep = 4000;

# source parameter
f0  = 7.0;
isz = 20; isx = floor(Int64, nx/2)+1;
zsrc = (isz-1)*deltaZ; xsrc = (isx-1)*deltaX;

# receivers
irz = [50, 50]; irx = [100, 200]

# ==============================================================================
#                        elastic model parameters
# ==============================================================================
# velocity model
vp = 3300.0; vs = vp / sqrt(3.0); density = 2800.0;
# arrays for elastic model parameter
lambda = zeros(nz+2, nx+2);
mu     = zeros(nz+2, nx+2);
rho    = zeros(nz+2, nx+2);

for j = 2 : nx+1
    for i = 2 : nz+1
           rho[i,j] = density
            mu[i,j] = density * vs * vs
        lambda[i,j] = density * (vp*vp - 2.0*vs*vs)
    end
end

# half grid elastic model parameters
hlambdaz    = 0.0
hmuz        = 0.0
hmux        = 0.0
hlambda2muz = 0.0
hrhozx      = 0.0

# ==============================================================================
#                       check stability conditions
# ==============================================================================
courantNumber = vp * dt * sqrt(1.0/deltaZ^2 + 1.0/deltaX^2)
if courantNumber > 1.0
   println("Courant number is: $courantNumber")
   error("time step is too large, simulation will be unstable")
end


# ==============================================================================
#                      prepare for CPML damping profile
# ==============================================================================
# thickness of CPML in number of grid points
nPointPml = 10;
thickPmlZ = nPointPml * deltaZ;
thickPmlX = nPointPml * deltaX;

# parameter for CPML profile
nPower   = 2.0
kmax     = 1.0
alphaMax = pi * f0

# damping profile, need to be improved  !!!
# (k, d, alpha, a, b) at grid point, # (hk, hd, halpha, ha, hb) at half grid point along Z direction
 kz=ones(nz);  dz=zeros(nz);  alphaz=zeros(nz);  az=zeros(nz);  bz=zeros(nz);
hkz=ones(nz); hdz=zeros(nz); halphaz=zeros(nz); haz=zeros(nz); hbz=zeros(nz);

# along X direction
 kx=zeros(nx);  dx=zeros(nx);  alphax=zeros(nx);  ax=zeros(nx);  bx=zeros(nx);
hkx=zeros(nx); hdx=zeros(nx); halphax=zeros(nx); hax=zeros(nx); hbx=zeros(nx);

# specify the obsorbing boundary, origin specify the computation domain, z
reflectCoef=1e-3;
d0Z = -(nPower + 1) * vp * log(reflectCoef) / (2.0 * thickPmlZ);
d0X = -(nPower + 1) * vp * log(reflectCoef) / (2.0 * thickPmlX);
locationZ=0.0; locationX=0.0; abscissa=0.0; normalizedAbscissa=0.0;
originTop = thickPmlZ; originBottom=(nz-1)*deltaZ - thickPmlZ;

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

    bz[i] = exp(- ( dz[i]/ kz[i] +  alphaz[i]) * dt)
   hbz[i] = exp(- (hdz[i]/hkz[i] + halphaz[i]) * dt)

    if (abs( dz[i]) > eps(Float32))
        az[i] =  dz[i] / ( kz[i]*( dz[i]+ kz[i]* alphaz[i])) * ( bz[i]-1.0)
    end
    if (abs(hdz[i]) > eps(Float32))
       haz[i] = hdz[i] / (hkz[i]*(hdz[i]+hkz[i]*halphaz[i])) * (hbz[i]-1.0)
    end
end

# damping in x direction
thickPmlX = nPointPml * deltaX; originLeft= thickPmlX; originRight =(nx-1)*deltaX - thickPmlX;
for j = 1 : nx
    locationX = deltaX * (j-1)
    # left damping profile at grid point
    abscissa = originLeft - locationX
    if (abscissa > 0.0)
       normalizedAbscissa = abscissa / thickPmlX
       dx[j] = d0X * normalizedAbscissa^nPower
       kx[j] = 1.0
       alphax[j] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end

    # left profile at half grid point
    abscissa = originLeft - (locationX + hdeltaX)
    if abscissa > 0.0
       normalizedAbscissa = abscissa / thickPmlX
       hdx[j] = d0X * normalizedAbscissa^npower
       hkx[j] = 1.0
       halphax[j] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end

    # right damping profile at grid point
    abscissa = locationX - originRight
    if abscissa > 0.0
       normalizedAbscissa = abscissa / thickPmlX
       dx[j] = d0X * normalizedAbscissa^nPower
       kx[j] = 1.0
       alphax[j] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end
    # bottom profile at half grid point
    abscissa = (locationX + hdeltaX) - originRight
    if abscissa > 0.0
       normalizedAbscissa = abscissa / thickPmlX
       hdx[j] = d0Z * normalizedAbscissa^npower
       hkx[j] = 1.0
       halphax[j] = alphaMax * (1.0 - normalizedAbscissa) + 0.1 * alphaMax
    end

    bx[j] = exp(- ( dx[j]/ kx[j] +  alphax[j]) * dt)
   hbx[j] = exp(- (hdx[j]/hkx[j] + halphax[j]) * dt)

    if (abs( dx[j]) > eps(Float32))
        ax[j] =  dx[j] / ( kx[j]*( dx[j]+ kx[j]* alphax[j])) * ( bx[j]-1.0)
    end
    if (abs(hdx[j]) > eps(Float32))
       hax[j] = hdx[j] / (hkx[j]*(hdx[j]+hkx[j]*halphax[j])) * (hbx[j]-1.0)
    end
end

# arrays for wave field variables
vz = zeros(nz+2, nx+2);
vx = zeros(nz+2, nx+2);
sigmazz = zeros(nz+2, nx+2);
sigmaxx = zeros(nz+2, nx+2);
sigmazx = zeros(nz+2, nx+2);

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
#                             start computation
# ==============================================================================







# crazy

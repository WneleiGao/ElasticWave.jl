using ElasticWave

nz = 173; nx = 267; # model size
ext= 20; # thickness of PML layers
iflag = 2; #when iflag=1, four side obsorbing boundary, iflag=2, free surface
dx = 10. ; dz = 10. ; dt = 1e-3; # grid size
tmax = 6; f0=10.; # recording time and dominant frequency of source wavelet
vp=3000.*ones(nz,nx);  vs=3000./sqrt(3)*ones(nz,nx); rho=2.5*ones(nz,nx);
# vp[150:end,:]=4000.;  vs[150:end,:]=4000./sqrt(3); #physical model
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

(lambda, mu) = vel2lm(vp, vs, rho);
spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1);
spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2);
wfd1 = InitWfd(nz, nx, ext, iflag, dt, 1);
wfd2 = InitWfd(nz, nx, ext, iflag, dt, 2);

# give random number to spt1
n = length(spt1.Vzz);
spt1.Vzz = randn(n);  spt1.Vzx = randn(n);  spt1.Vxz = randn(n);  spt1.Vxx = randn(n);
spt1.Tzzz = randn(n); spt1.Tzzx = randn(n); spt1.Txxz = randn(n); spt1.Txxx = randn(n);
spt1.Txzz = randn(n); spt1.Txzx = randn(n);
OneStepForward!(spt2, spt1, fidMtx);

#convert spt to wfd
if iflag == 1
   Nz = nz + 2*ext
   zu = ext
elseif iflag == 2
   Nz = nz +   ext
   zu = 0
end
Nx = nx + 2*ext
tmp = reshape(spt1.Vzz+spt1.Vzx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd1.Vz = vec(tmp);
tmp = reshape(spt1.Vxz+spt1.Vxx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd1.Vx = vec(tmp);
tmp = reshape(spt1.Tzzz+spt1.Tzzx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd1.Tzz = vec(tmp);
tmp = reshape(spt1.Txxz+spt1.Txxx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd1.Txx = vec(tmp);
tmp = reshape(spt1.Txzz+spt1.Txzx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd1.Txz = vec(tmp);

# save the boundary of spt1
lb = (4*nz+(nx-4)*4)*5;
bd1 = zeros(lb);
tmp=reshape(wfd1.Vz,nz,nx);
il=   1; iu=   2*nz    ; bd1[il:iu]=vec(tmp[:,1:2]);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,end-1:end]);
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[1:2,3:nx-2])');
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[end-1:end,3:nx-2])');

tmp=reshape(wfd1.Vx,nz,nx);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,1:2]);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,end-1:end]);
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[1:2,3:nx-2])');
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[end-1:end,3:nx-2])');

tmp=reshape(wfd1.Tzz,nz,nx);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,1:2]);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,end-1:end]);
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[1:2,3:nx-2])');
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[end-1:end,3:nx-2])');

tmp=reshape(wfd1.Txx,nz,nx);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,1:2]);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,end-1:end]);
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[1:2,3:nx-2])');
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[end-1:end,3:nx-2])');

tmp=reshape(wfd1.Txz,nz,nx);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,1:2]);
il=iu+1; iu=iu+2*nz    ; bd1[il:iu]=vec(tmp[:,end-1:end]);
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[1:2,3:nx-2])');
il=iu+1; iu=iu+2*(nx-4); bd1[il:iu]=vec((tmp[end-1:end,3:nx-2])');

bv = zeros(bd1);
extractBoundary!(bv, spt1);


tmp = reshape(spt2.Vzz+spt2.Vzx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd2.Vz = vec(tmp);
tmp = reshape(spt2.Vxz+spt2.Vxx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd2.Vx = vec(tmp);
tmp = reshape(spt2.Tzzz+spt2.Tzzx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd2.Tzz = vec(tmp);
tmp = reshape(spt2.Txxz+spt2.Txxx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd2.Txx = vec(tmp);
tmp = reshape(spt2.Txzz+spt2.Txzx, Nz, Nx);
tmp = tmp[zu+1:zu+nz, ext+1:ext+nx]; wfd2.Txz = vec(tmp);

vz2 = reshape(wfd2.Vz, nz, nx);   vz1 = zeros(vz2);
vx2 = reshape(wfd2.Vx, nz, nx);   vx1 = zeros(vx2);
Tzz2 = reshape(wfd2.Tzz, nz, nx); Tzz1 = zeros(Tzz2);
Txx2 = reshape(wfd2.Txx, nz, nx); Txx1 = zeros(Txx2);
Txz2 = reshape(wfd2.Txz, nz, nx); Txz1 = zeros(Txz2);

vz = reshape(wfd1.Vz, nz, nx);
vx = reshape(wfd1.Vx, nz, nx);
Tzz = reshape(wfd1.Tzz, nz, nx);
Txx = reshape(wfd1.Txx, nz, nx);
Txz = reshape(wfd1.Txz, nz, nx);

a1 =  9/8; a2 = -1/24;
for j = 3:nx-2
    for i = 3 : nz-2
        c  = dt*(mu[i,j]+mu[i+1,j]+mu[i,j+1]+mu[i+1,j+1])/4
        c1 = c*a1/dz; c2 = c*a2/dz;
        t1 = c1*(vx2[i+1,j]-vx2[i  ,j])
        t2 = c2*(vx2[i+2,j]-vx2[i-1,j])
        c1 = c*a1/dx; c2 = c*a2/dx;
        t3 = c1*(vz2[i,j+1]-vz2[i,j  ])
        t4 = c2*(vz2[i,j+2]-vz2[i,j-1])
        Txz1[i,j] = Txz2[i,j]-(t1+t2+t3+t4)
    end
end

il = (4*nz+(nx-4)*4)*4;
Txz1[:,1:2] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
Txz1[:,end-1:end] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
Txz1[1:2,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))'; il = il+2*(nx-4);
Txz1[end-1:end,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))';

for j = 3:nx-2
    for i = 3 : nz-2
        c  = dt*(lambda[i,j]+2*mu[i,j])
        c1 = c*a1/dx; c2 = c*a2/dx;
        t1 = c1*(vx2[i,j  ]-vx2[i,j-1])
        t2 = c2*(vx2[i,j+1]-vx2[i,j-2])
        c  = dt*lambda[i,j]
        c1 = c*a1/dz; c2 = c*a2/dz;
        t3 = c1*(vz2[i  ,j]-vz2[i-1,j])
        t4 = c2*(vz2[i+1,j]-vz2[i-2,j])
        Txx1[i,j] = Txx2[i,j]-(t1+t2+t3+t4)
    end
end
il = (4*nz+(nx-4)*4)*3;
Txx1[:,1:2] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
Txx1[:,end-1:end] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
Txx1[1:2,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))'; il = il+2*(nx-4);
Txx1[end-1:end,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))';


for j = 3:nx-2
    for i = 3 : nz-2
        c  = dt*(lambda[i,j]+2*mu[i,j])
        c1 = c*a1/dz; c2 = c*a2/dz;
        t1 = c1*(vz2[i  ,j]-vz2[i-1,j])
        t2 = c2*(vz2[i+1,j]-vz2[i-2,j])
        c  = dt*lambda[i,j]
        c1 = c*a1/dx; c2 = c*a2/dx;
        t3 = c1*(vx2[i,j  ]-vx2[i,j-1])
        t4 = c2*(vx2[i,j+1]-vx2[i,j-2])
        Tzz1[i,j] = Tzz2[i,j]-(t1+t2+t3+t4)
    end
end
il = (4*nz+(nx-4)*4)*2;
Tzz1[:,1:2] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
Tzz1[:,end-1:end] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
Tzz1[1:2,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))'; il = il+2*(nx-4);
Tzz1[end-1:end,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))';

# Vx
for j = 3:nx-2
    for i = 3 : nz-2
        c  = 2*dt/(rho[i,j]+rho[i,j+1])
        c1 = c*a1/dx; c2 = c*a2/dx;
        t1 = c1*(Txx1[i,j+1]-Txx1[i,j  ])
        t2 = c2*(Txx1[i,j+2]-Txx1[i,j-1])
        c1 = c*a1/dz; c2 = c*a2/dz;
        t3 = c1*(Txz1[i  ,j]-Txz1[i-1,j])
        t4 = c2*(Txz1[i+1,j]-Txz1[i-2,j])
        vx1[i,j] = vx2[i,j]-(t1+t2+t3+t4)
    end
end
il = (4*nz+(nx-4)*4);
vx1[:,1:2] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
vx1[:,end-1:end] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
vx1[1:2,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))'; il = il+2*(nx-4);
vx1[end-1:end,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))';

# vz
for j = 3:nx-2
    for i = 3 : nz-2
        c  = 2*dt/(rho[i+1,j]+rho[i,j])
        c1 = c*a1/dz; c2 = c*a2/dz;
        t1 = c1*(Tzz1[i+1,j]-Tzz1[i  ,j])
        t2 = c2*(Tzz1[i+2,j]-Tzz1[i-1,j])
        c1 = c*a1/dx; c2 = c*a2/dx;
        t3 = c1*(Txz1[i,j  ]-Txz1[i,j-1])
        t4 = c2*(Txz1[i,j+1]-Txz1[i,j-2])
        vz1[i,j] = vz2[i,j]-(t1+t2+t3+t4)
    end
end
il = 0;
vz1[:,1:2] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
vz1[:,end-1:end] = reshape(bd1[il+1:il+2*nz], nz, 2); il = il+2*nz;
vz1[1:2,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))'; il = il+2*(nx-4);
vz1[end-1:end,3:end-2] = (reshape(bd1[il+1:il+2*(nx-4)],(nx-4),2))';


# =====test write boundary to a binary file =====
using ElasticWave
nz = 100; nx = 100; # model size
ext= 20; # thickness of PML layers
iflag = 2; #when iflag=1, four side obsorbing boundary, iflag=2, free surface
dx = 10. ; dz = 10. ; dt = 1e-3; # grid size
tmax = 3.; f0=10.; # recording time and dominant frequency of source wavelet
vp=3000.*ones(nz,nx);  vs=3000./sqrt(3)*ones(nz,nx); rho=2.5e-4*ones(nz,nx);
vp[150:end,:]=4000.;  vs[150:end,:]=4000./sqrt(3);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

isx = 50; isz = 1; ot=0.
flags = vec([false false true true false]); # add source to normal stress, simulate explosive sources
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

root = homedir();
path = joinpath(root, "Desktop/bv.bin")
MultiStepForward(path, src, fidMtx, tmax)
bv = ReadWfdBoundary(path)

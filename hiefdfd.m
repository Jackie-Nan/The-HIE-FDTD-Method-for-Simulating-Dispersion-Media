
tic
light=3*1e+8;
e0=1/(36*pi)*1e-9;
u0=4*pi*1e-7;
wavelength=0.15;
%%
px1=1;
px2=11;
obx1=20;  
obx2=39; 
px3=50;
px4=60;
%%
py1=1;
py2=11;
oby1=20;  
oby2=39;
py3=50;
py4=60;
%%
pz1=1;
pz2=11;
obz1=20;  
subz11=20;
subz12=22;
debyz11=23;
debyz12=24;
durdz11=25;
durdz12=26;
Lorentz11=27;
Lorentz12=28;
subz21=29;
subz22=31;
debyz21=32;
debyz22=33;
durdz21=34;
durdz22=35;
Lorentz21=36;
Lorentz22=37;
subz31=38;
subz32=40;
debyz31=41;
debyz32=42;
durdz31=43;
durdz32=44;
Lorentz31=45;
Lorentz32=46;
subz41=47;
subz42=49;
obz2=49;
pz3=60;
pz4=70;
%%%
zdeta=zeros(pz4,1)+wavelength/10;
for p=subz11:1:subz12
    zdeta(p)=0.015;
end
for p=debyz11:1: debyz12
    zdeta(p)=0.0005;
end
for p=durdz11:1: durdz12
    zdeta(p)=0.0005;
end
for p=Lorentz11:1: Lorentz12
    zdeta(p)=0.0005;
end
for p=subz21:1:subz22
    zdeta(p)=0.015;
end
for p=debyz21:1: debyz22
    zdeta(p)=0.0005;
end
for p=durdz21:1: durdz22
    zdeta(p)=0.0005;
end
for p=Lorentz21:1: Lorentz22
    zdeta(p)=0.0005;
end
for p=subz31:1:subz32
    zdeta(p)=0.015;
end
for p=debyz31:1: debyz32
    zdeta(p)=0.0005;
end
for p=durdz31:1: durdz32
    zdeta(p)=0.0005;
end
for p=Lorentz31:1: Lorentz32
    zdeta(p)=0.0005;
end
for p=subz41:1:subz42
    zdeta(p)=0.015;
end

xdeta=zeros(px4,1)+wavelength/10;
ydeta=zeros(py4,1)+wavelength/10;
%%
zde=min(zdeta);
xde=min(xdeta);
yde=min(ydeta);
%%
tdeta=1/sqrt((1/xde)^2+(1/yde)^2)/light;
stepnumber=300; 
%&
omiga=2*pi*light/wavelength;  
tt=1*1e-9;  
t0=1*1e-9;
%
i1=15;
i2=45;
j1=15;
j2=45;
k1=15;
k2=55;
%
fieldx=30;
fieldy=30;
fieldz=35;
fieldx1=30;
fieldy1=18;
fieldz1=35;
fieldx2=30;
fieldy2=43;
fieldz2=35;

%
n01=2.40483;
ncfsl=px2-px1;
ncfsr=px4-px3;
ncfsf=py2-py1;
ncfsb=py4-py3;
ncfsd=pz2-pz1;
ncfsu=pz4-pz3;
detamaxl=2.12;
detamaxr=2.12;
detamaxf=2.12;
detamaxb=2.12;
detamaxd=2.12;
detamaxu=2.12;
kmaxl=15;
kmaxr=15;
kmaxf=15;
kmaxb=15;
kmaxd=15;
kmaxu=15;
mmil=4;
mmir=4;
mmif=4;
mmib=4;
mmid=4;
mmiu=4;
afa=0.08;
%%%%%%%%%
ax=zeros(px4+1,1);
kx=ax+1;
axm=zeros(px4,1);
kxm=axm+1;
bx=zeros(px4+1,1);
bxm=zeros(px4,1);
for m=px1:1:px2-1
    hao=detamaxl*(abs(m-px2)/ncfsl)^mmil;
    kx(m)=1+(kmaxl-1)*(abs(m-px2)/ncfsl)^mmil;
    bx(m)=exp(-(hao/kx(m)+afa)*tdeta/e0);
    ax(m)=hao/(hao*kx(m)+kx(m)^2*afa)*(bx(m)-1);
     hao=detamaxl*(abs(m-px2+1/2)/ncfsl)^mmil;
     kxm(m)=1+(kmaxl-1)*(abs(m-px2+1/2)/ncfsl)^mmil;
     bxm(m)=exp(-(hao/kxm(m)+afa)*tdeta/e0);
    axm(m)=hao/(hao*kxm(m)+kxm(m)^2*afa)*(bxm(m)-1);
end
for m=px3+2:1:px4+1
      hao=detamaxr*(abs(m-px3-1)/ncfsr)^mmir;
    kx(m)=1+(kmaxr-1)*(abs(m-px3-1)/ncfsr)^mmir;
    bx(m)=exp(-(hao/kx(m)+afa)*tdeta/e0);
    ax(m)=hao/(hao*kx(m)+kx(m)^2*afa)*(bx(m)-1);
     hao=detamaxr*(abs(m-px3-1-1/2)/ncfsr)^mmir;
     kxm(m-1)=1+(kmaxr-1)*(abs(m-px3-1-1/2)/ncfsr)^mmir;
     bxm(m-1)=exp(-(hao/kxm(m-1)+afa)*tdeta/e0);
    axm(m-1)=hao/(hao*kxm(m-1)+kxm(m-1)^2*afa)*(bxm(m-1)-1);
    
end
%%%%%%%%%%%%
ay=zeros(py4+1,1);
ky=ay+1;
aym=zeros(py4,1);
kym=aym+1;
by=zeros(py4+1,1);
bym=zeros(py4,1);
for n=py1:1:py2-1    
      hao=detamaxf*(abs(n-py2)/ncfsf)^mmif;
    ky(n)=1+(kmaxf-1)*(abs(n-py2)/ncfsf)^mmif;
    by(n)=exp(-(hao/ky(n)+afa)*tdeta/e0);
    ay(n)=hao/(hao*ky(n)+ky(n)^2*afa)*(by(n)-1);
     hao=detamaxf*(abs(n-py2+1/2)/ncfsf)^mmif;
     kym(n)=1+(kmaxf-1)*(abs(n-py2+1/2)/ncfsf)^mmif;
     bym(n)=exp(-(hao/kym(n)+afa)*tdeta/e0);
    aym(n)=hao/(hao*kym(n)+kym(n)^2*afa)*(bym(n)-1);
end
for n=py3+2:1:py4+1
    hao=detamaxb*(abs(n-py3-1)/ncfsb)^mmib;
    ky(n)=1+(kmaxb-1)*(abs(n-py3-1)/ncfsb)^mmib;
    by(n)=exp(-(hao/ky(n)+afa)*tdeta/e0);
    ay(n)=hao/(hao*ky(n)+ky(n)^2*afa)*(by(n)-1);
     hao=detamaxb*(abs(n-py3-1-1/2)/ncfsb)^mmib;
     kym(n-1)=1+(kmaxb-1)*(abs(n-py3-1-1/2)/ncfsb)^mmib;
     bym(n-1)=exp(-(hao/kym(n-1)+afa)*tdeta/e0);
    aym(n-1)=hao/(hao*kym(n-1)+kym(n-1)^2*afa)*(bym(n-1)-1);
end
%%%%%%%%%%
az=zeros(pz4+1,1);
kz=az+1;
azm=zeros(pz4,1);
kzm=azm+1;
bz=zeros(pz4+1,1);
bzm=zeros(pz4,1);
for p=pz1:1:pz2-1
    
     hao=detamaxd*(abs(p-pz2)/ncfsd)^mmid;
    kz(p)=1+(kmaxd-1)*(abs(p-pz2)/ncfsd)^mmid;
    bz(p)=exp(-(hao/kz(p)+afa)*tdeta/e0);
    az(p)=hao/(hao*kz(p)+kz(p)^2*afa)*(bz(p)-1);
     hao=detamaxd*(abs(p-pz2+1/2)/ncfsd)^mmid;
     kzm(p)=1+(kmaxd-1)*(abs(p-pz2+1/2)/ncfsd)^mmid;
     bzm(p)=exp(-(hao/kzm(p)+afa)*tdeta/e0);
    azm(p)=hao/(hao*kzm(p)+kzm(p)^2*afa)*(bzm(p)-1);
end
for p=pz3+2:1:pz4+1
      hao=detamaxu*(abs(p-pz3-1)/ncfsu)^mmiu;
    kz(p)=1+(kmaxu-1)*(abs(p-pz3-1)/ncfsu)^mmiu;
    bz(p)=exp(-(hao/kz(p)+afa)*tdeta/e0);
    az(p)=hao/(hao*kz(p)+kz(p)^2*afa)*(bz(p)-1);
     hao=detamaxu*(abs(p-pz3-1-1/2)/ncfsu)^mmiu;
     kzm(p-1)=1+(kmaxu-1)*(abs(p-pz3-1-1/2)/ncfsu)^mmiu;
     bzm(p-1)=exp(-(hao/kzm(p-1)+afa)*tdeta/e0);
    azm(p-1)=hao/(hao*kzm(p-1)+kzm(p-1)^2*afa)*(bzm(p-1)-1);
end

%
cp=1;
cq=tdeta/u0;
coff= tdeta^2/(4*e0*u0);
%
segma= zeros(pz4+1,1);
%
for p=subz11:1:subz12
segma(p)=0.002;
end
for p=subz21:1:subz22
segma(p)=0.002;
end
for p=subz31:1:subz32
segma(p)=0.002;
end
for p=subz41:1:subz42
segma(p)=0.002;
end
%
er= zeros(pz4+1,1)+1;
for p=subz11:1:subz12
er(p)=2.45;
end
for p=subz21:1:subz22
er(p)=2.45;
end
for p=subz31:1:subz32
er(p)=2.45;
end
for p=subz41:1:subz42
er(p)=2.45;
end
%
eeff=e0*er;
%
coffe=zeros(pz4+1,1);
for p=1:1:pz4+1
coffe(p)= tdeta^2/(4*e0*er(p,1)*u0);
end

%
kp2=zeros(pz4+1,1);  
kp3=zeros(pz4+1,1);  
De7=zeros(pz4+1,1);  
De8=zeros(pz4+1,1);  
for p=1:1: pz4+1
De7(p,1)=1+segma(p)/2*tdeta/(e0*er(p,1));
De8(p,1)=1- segma(p)/2*tdeta/(e0*er(p,1));
end 
De9=zeros(pz4+1,1);  

De10=zeros(pz4+1,1);  
De11=zeros(pz4+1,1);  

betaL1= zeros(pz4+1,1);  
betaL2= zeros(pz4+1,1);  
betap1=zeros(pz4+1,1);  

%%%
omiga0=2*pi*10^9;  
vc=4*10^9;  
es=zeros(pz4+1,1)+1;
einfi=zeros(pz4+1,1)+1;
for p=Lorentz11:1: Lorentz12  
es(p,1)=2;
einfi(p,1)= 1;
eeff(p,1)=e0*einfi(p,1);

end
for p=Lorentz21:1: Lorentz22  
es(p,1)=2;
einfi(p,1)= 1;
eeff(p,1)=e0*einfi(p,1);
end
for p=Lorentz31:1: Lorentz32  
es(p,1)=2;
einfi(p,1)= 1;
eeff(p,1)=e0*einfi(p,1);
end
edeta=es-einfi; 
%
%s
for p= Lorentz11:1: Lorentz12
coffe(p)= tdeta^2/(4* eeff(p,1)*u0);
end
for p= Lorentz21:1: Lorentz22
coffe(p)= tdeta^2/(4* eeff(p,1)*u0);
end
for p= Lorentz31:1: Lorentz32
coffe(p)= tdeta^2/(4* eeff(p,1)*u0);
end

for p= Lorentz11:1: Lorentz12
kp2(p)=(2-omiga0^2*tdeta^2)/(1+vc*tdeta);
kp3(p)=( vc*tdeta -1)/(1+vc*tdeta);
end
for p= Lorentz21:1: Lorentz22
kp2(p)=(2-omiga0^2*tdeta^2)/(1+vc*tdeta);
kp3(p)=( vc*tdeta -1)/(1+vc*tdeta);
end

for p= Lorentz31:1: Lorentz32
kp2(p)=(2-omiga0^2*tdeta^2)/(1+vc*tdeta);
kp3(p)=( vc*tdeta -1)/(1+vc*tdeta);
end

betap2=zeros(pz4+1,1);  %
for p= Lorentz11:1: Lorentz12
betap2(p,1)=e0*edeta(p,1)*omiga0^2*tdeta^2/(1+vc*tdeta);
De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+betap2(p,1)/4/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1);
De9(p,1)=-(kp2(p)+1)/2*tdeta/eeff(p,1);
De10(p,1)= betap2(p,1)/4/eeff(p,1);
De11(p,1)=tdeta/eeff(p,1) *kp3(p,1)/2;
betaL1(p,1)=betap2(p,1)/2/tdeta;
betaL2(p,1)=betap2(p,1)/2/tdeta;
end
for p= Lorentz21:1: Lorentz22
betap2(p,1)=e0*edeta(p,1)*omiga0^2*tdeta^2/(1+vc*tdeta);
De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+betap2(p,1)/4/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1);
De9(p,1)=-(kp2(p)+1)/2*tdeta/eeff(p,1);
De10(p,1)= betap2(p,1)/4/eeff(p,1);
De11(p,1)=tdeta/eeff(p,1) *kp3(p,1)/2;
betaL1(p,1)=betap2(p,1)/2/tdeta;
betaL2(p,1)=betap2(p,1)/2/tdeta;
end

for p= Lorentz31:1: Lorentz32
betap2(p,1)=e0*edeta(p,1)*omiga0^2*tdeta^2/(1+vc*tdeta);
De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+betap2(p,1)/4/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1);
De9(p,1)=-(kp2(p)+1)/2*tdeta/eeff(p,1);
De10(p,1)= betap2(p,1)/4/eeff(p,1);
De11(p,1)=tdeta/eeff(p,1) *kp3(p,1)/2;
betaL1(p,1)=betap2(p,1)/2/tdeta;
betaL2(p,1)=betap2(p,1)/2/tdeta;
end

%
omigap=2*pi*10^9; 
vc=4*10^9; 
for p=durdz11:1: durdz12
einfi(p,1)=1;
eeff(p,1)=e0*einfi(p,1);
coffe(p)= tdeta^2/(4*eeff(p,1)*u0);
end
for p=durdz21:1: durdz22
einfi(p,1)=1;
eeff(p,1)=e0*einfi(p,1);
coffe(p)= tdeta^2/(4*eeff(p,1)*u0);
end
for p=durdz31:1: durdz32
einfi(p,1)=1;
eeff(p,1)=e0*einfi(p,1);
coffe(p)= tdeta^2/(4*eeff(p,1)*u0);
end


for p= durdz11:1: durdz12
kp2(p)= (2-vc*tdeta)/(2+vc*tdeta);
end
for p= durdz21:1: durdz22
kp2(p)= (2-vc*tdeta)/(2+vc*tdeta);
end
for p= durdz31:1: durdz32
kp2(p)= (2-vc*tdeta)/(2+vc*tdeta);
end

for p= durdz11:1: durdz12
De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+ e0*omigap^2*tdeta/(2+vc*tdeta)/2*tdeta/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1)- e0*omigap^2*tdeta/(2+vc*tdeta)/2*tdeta/eeff(p,1);
De9(p,1)=-( (2-vc*tdeta)/(2+vc*tdeta)+1)/2*tdeta/eeff(p,1);

betaL1(p,1)= e0*omigap^2*tdeta/(2+vc*tdeta) ;
betap1(p,1)= e0*omigap^2*tdeta/(2+vc*tdeta)
end
for p= durdz21:1: durdz22
De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+ e0*omigap^2*tdeta/(2+vc*tdeta)/2*tdeta/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1)- e0*omigap^2*tdeta/(2+vc*tdeta)/2*tdeta/eeff(p,1);
De9(p,1)=-( (2-vc*tdeta)/(2+vc*tdeta)+1)/2*tdeta/eeff(p,1);

betaL1(p,1)= e0*omigap^2*tdeta/(2+vc*tdeta) ;
betap1(p,1)= e0*omigap^2*tdeta/(2+vc*tdeta)
end
for p= durdz31:1: durdz32
De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+ e0*omigap^2*tdeta/(2+vc*tdeta)/2*tdeta/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1)- e0*omigap^2*tdeta/(2+vc*tdeta)/2*tdeta/eeff(p,1);
De9(p,1)=-( (2-vc*tdeta)/(2+vc*tdeta)+1)/2*tdeta/eeff(p,1);

betaL1(p,1)= e0*omigap^2*tdeta/(2+vc*tdeta) ;
betap1(p,1)= e0*omigap^2*tdeta/(2+vc*tdeta);
end

%
vc=1*10^8; 
tao=1/vc;
es=zeros(pz4+1,1)+1;
for p=debyz11:1: debyz12 
es(p,1)=2;
einfi(p,1)=1;
eeff(p,1)=e0*einfi(p,1);
coffe(p)= tdeta^2/(4*eeff(p,1)*u0);
end
for p=debyz21:1: debyz22  
es(p,1)=2;
einfi(p,1)= 1;
eeff(p,1)=e0*einfi(p,1);
coffe(p)= tdeta^2/(4*eeff(p,1)*u0);
end
for p=debyz31:1: debyz32  
es(p,1)=2;
einfi(p,1)= 1;
eeff(p,1)=e0*einfi(p,1);
coffe(p)= tdeta^2/(4*eeff(p,1)*u0);
end

edeta=es-einfi;


for p= debyz11:1: debyz12
kp2(p)= (2*tao- tdeta)/(2*tao+tdeta);
end
for p= debyz21:1: debyz22
kp2(p)= (2*tao- tdeta)/(2*tao+tdeta);
end
for p= debyz31:1: debyz32
kp2(p)= (2*tao- tdeta)/(2*tao+tdeta);
end
for p= debyz11:1: debyz12

De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+ (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/2/tdeta*tdeta/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1)- (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/eeff(p,1);
De9(p,1)=-(kp2(p,1)+1)/2*tdeta/eeff(p,1);

betaL1(p,1)= (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/tdeta;
betap1(p,1)= (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/tdeta;
end
for p= debyz21:1: debyz22

De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+ (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/2/tdeta*tdeta/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1)- (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/eeff(p,1);
De9(p,1)=-(kp2(p,1)+1)/2*tdeta/eeff(p,1);

betaL1(p,1)= (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/tdeta;
betap1(p,1)= (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/tdeta;
end
for p= debyz31:1: debyz32

De7(p,1)=1+segma(p)/2*tdeta/eeff(p,1)+ (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/2/tdeta*tdeta/eeff(p,1);
De8(p,1)=1- segma(p)/2*tdeta/eeff(p,1)- (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/eeff(p,1);
De9(p,1)=-(kp2(p,1)+1)/2*tdeta/eeff(p,1);

betaL1(p,1)= (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/tdeta;
betap1(p,1)= (2*e0*edeta(p,1)* tdeta/(2*tao+tdeta))/tdeta;
end

%%
matfree=zeros(pz4-1,pz4-1);
for m=1:1:pz4-1
  matfree(m,m)=1+coff/((zdeta(m+1)+zdeta(m))/2)/zdeta(m+1)*(1/kz(m+1)+az(m+1))*(1/kzm(m+1)+azm(m+1))+coff/((zdeta(m+1)+zdeta(m))/2)/zdeta(m)*(1/kz(m+1)+az(m+1))*(1/kzm(m)+azm(m));
 end
 for m=2:1:pz4-1
     matfree(m,m-1)=-coff/((zdeta(m+1)+zdeta(m))/2)/zdeta(m)*(1/kz(m+1)+az(m+1))*(1/kzm(m)+azm(m));
 end
 for m=1:1:pz4-2
     matfree(m,m+1)=-coff/((zdeta(m+1)+zdeta(m))/2)/zdeta(m+1)*(1/kz(m+1)+az(m+1))*(1/kzm(m+1)+azm(m+1));
 end
matfree=inv(matfree);

%
mat=zeros(pz4-1,pz4-1);

for m= 1:1:pz4-1
mat(m,m)=De7(m,1)+coffe(p,1)/((zdeta(m+1)+zdeta(m))/2)/zdeta(m+1)*(1/kz(m+1)+az(m+1))*(1/kzm(m+1)+azm(m+1))+coffe(p,1)/((zdeta(m+1)+zdeta(m))/2)/zdeta(m)*(1/kz(m+1)+az(m+1))*(1/kzm(m)+azm(m));
 end

 for m=2:1:pz4-1
     mat(m,m-1)=-coffe(p,1)/((zdeta(m+1)+zdeta(m))/2)/zdeta(m)*(1/kz(m+1)+az(m+1))*(1/kzm(m)+azm(m));
 end

for m=1:1:pz4-2
     mat(m,m+1)=-coffe(p,1)/((zdeta(m+1)+zdeta(m))/2)/zdeta(m+1)*(1/kz(m+1)+az(m+1))*(1/kzm(m+1)+azm(m+1));
 end
mat=inv(mat);
so1=zeros(2,1);  
so2=zeros(2,1); 
%
ex1=zeros(px4,py4+1,pz4+1);
ey1=zeros(px4+1,py4,pz4+1);
ez1=zeros(px4+1,py4+1,pz4);
ex0=zeros(px4,py4+1,pz4+1);
ey0=zeros(px4+1,py4,pz4+1);
ez0=zeros(px4+1,py4+1,pz4);
Jpx1=zeros(px4,py4+1,pz4+1);
Jpy1=zeros(px4+1,py4,pz4+1);
Jpz1=zeros(px4+1,py4+1,pz4);
Jpx0=zeros(px4,py4+1,pz4+1);
Jpy0=zeros(px4+1,py4,pz4+1);
Jpz0=zeros(px4+1,py4+1,pz4);

hz1=zeros(px4,py4,pz4+1);
hy1=zeros(px4,py4+1,pz4);
hx1=zeros(px4+1,py4,pz4);
ex2=ex1;
ey2=ey1;
ez2=ez1;
hx2=hx1;
hy2=hy1;
hz2=hz1;
%
ezi=ez1;
ezi1=ez1;
hxi=hx1;

mxy1=hx1;
mxz1=hx1;
mxy2=mxy1;
mxz2=mxz1;
myz1=hy1;
myx1=hy1;
myz2=myz1;
myx2=myx1;
mzy1=hz1;
mzx1=hz1;
mzy2=mzy1;
mzx2=mzx1;
pxy1=ex1;
pxz1=ex1;
pxy2=pxy1;
pxz2=pxz1;
pyx1=ey1;
pyz1=ey1;
pyx2=pyx1;
pyz2=pyz1;
pzx1=ez1;
pzy1=ez1;
pzx2=pzx1;
pzy2=pzy1;

incidentstart=j1-3;
incidentend=j2+3;
ezin1=zeros(incidentend+1,1);
ezin2=zeros(incidentend+1,1);
hxin1=zeros(incidentend,1);
hxin2=zeros(incidentend,1);

resulthex=zeros(1,1);
resulthey=zeros(1,1);
resulthez=zeros(1,1);
resulthhx=zeros(1,1);
resulthhy=zeros(1,1);
resulthhz=zeros(1,1);
resulthex1=zeros(1,1);
resulthey1=zeros(1,1);
resulthez1=zeros(1,1);
resulthhx1=zeros(1,1);
resulthhy1=zeros(1,1);
resulthhz1=zeros(1,1);
resulthex2=zeros(1,1);
resulthey2=zeros(1,1);
resulthez2=zeros(1,1);
resulthhx2=zeros(1,1);
resulthhy2=zeros(1,1);
resulthhz2=zeros(1,1);


source=zeros(1,1);
xzh=zeros(1,1); 
%
for t=1:1:stepnumber
    t
    xzh(t)=t*tdeta*1e+9; 
% 
  yv=exp(-4*pi*(t*tdeta-t0)^2/tt^2);
    for i=incidentstart+1:1:incidentend
    ezin2(i)=ezin1(i)-tdeta/e0/((ydeta(i)+ydeta(i-1))/2)*(hxin1(i)-hxin1(i-1));
  end
   ezin2(incidentend+1)=ezin1(incidentend)+(light*tdeta-ydeta(incidentend))/(light*tdeta+ydeta(incidentend))*(ezin2(incidentend)-ezin1(incidentend+1));
   ezin2(incidentstart)=yv;
  for i=incidentstart:1:incidentend
    hxin2(i)=hxin1(i)-tdeta/u0/ydeta(i)*(ezin2(i+1)-ezin2(i));
  end
  source(t)=ezin2(incidentstart+5);  
 
 for m=2:1:px4
for n=2:1:py4
if m<obx1|m>obx2|n<oby1|n>oby2 
for p=1:1:pz4
            ez2(m,n,p)= ez1(m,n,p)+tdeta/e0*((hy1(m,n,p)-hy1(m-1,n,p))/ ((xdeta(m)+xdeta(m-1))/2)/kx(m)-(hx1(m,n,p)-hx1(m,n-1,p))/ky(n)/ ((ydeta(n)+ydeta(n-1))/2))+ tdeta/ e0 *(pzx1(m,n,p)-pzy1(m,n,p));
end
           else 
for p=1:1:pz4
            ez2(m,n,p)=De8(p,1)*ez1(m,n,p)+tdeta/eeff(p,1)*((hy1(m,n,p)-hy1(m-1,n,p))/ ((xdeta(m)+xdeta(m-1))/2)/kx(m)-(hx1(m,n,p)-hx1(m,n-1,p))/ky(n)/ ((ydeta(n)+ydeta(n-1))/2))+ tdeta/ eeff(p,1) *(pzx1(m,n,p)-pzy1(m,n,p));
ez2(m,n,p)=ez2(m,n,p)+De9(p,1)*Jpz1(m,n,p)+ De10(p,1) * ez0(m,n,p)-De11(p,1) * Jpz0(m,n,p);
ez2(m,n,p)=1/De7(p,1)*ez2(m,n,p);
         end 
  end
end
     end
%
  for m=i1:1:i2+1
     for p=k1:1:k2
         ez2(m,j2+1,p)= ez1(m,j2+1,p)+tdeta/ e0 *((hy1(m,j2+1,p)-hy1(m-1,j2+1,p))/ ((xdeta(m)+xdeta(m-1))/2)-(hx1(m,j2+1,p)-hx1(m,j2+1-1,p))/((ydeta(j2+1)+ydeta(j2))/2))- tdeta/e0 /((ydeta(j2+1)+ydeta(j2))/2)*hxin1(j2+1);

         ez2(m,j1,p)= ez1(m,j1,p)+tdeta/ e0 *((hy1(m,j1,p)-hy1(m-1,j1,p))/ ((xdeta(m)+xdeta(m-1))/2)-(hx1(m,j1,p)-hx1(m,j1-1,p))/ ((ydeta(j1)+ydeta(j1-1))/2))+ tdeta/e0/((ydeta(j1)+ydeta(j1-1))/2)*hxin1(j1-1);
     end
 end


%
for m=obx1:1:obx2
    for n=oby1:1:oby2
        for p=1:1:pz4
Jpz2(m,n,p)=kp2(p,1)*Jpz1(m,n,p)+kp3(p,1)*Jpz0(m,n,p)+betaL1(p,1)* ez2(m,n,p)- betaL2(p,1)* ez0(m,n,p) +betap1(p,1)*ez1(m,n,p);
end
end
end
 %
for m=px1:1:px4
    for n=py1:1:py4+1
        for p=pz1:1:pz4
            myx2(m,n,p)=bxm(m)*myx1(m,n,p)+axm(m)/xdeta(m)*(ez2(m+1,n,p)-ez2(m,n,p));
        end
    end
end
    for m=px1:1:px4+1
    for n=py1:1:py4
        for p=pz1:1:pz4
           mxy2(m,n,p)=bym(n)*mxy1(m,n,p)+aym(n)/ydeta(n)*(ez2(m,n+1,p)-ez2(m,n,p));
        end
    end
    end
% 
%
for n=j1:1:j2+1
    for p=k1:1:k2
       ezi(i1-1,n,p)=-cq*ezin2(n)/xdeta(i1-1);
       ezi(i2+1,n,p)=cq*ezin2(n)/xdeta(i2+1);
    end
end 

   for m=1:1:px4
        for n=2:1:py4 
if m<obx1|m>obx2|n<oby1|n>oby2

       for p=2:1:pz4
            co1=tdeta/e0/((ydeta(n)+ydeta(n-1))/2)*(1/ky(n))*(hz1(m,n,p)-hz1(m,n-1,p))- tdeta^2/(2*u0*e0*xdeta(m))/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(1/kxm(m)+axm(m))*(ez2(m+1,n,p)-ez2(m,n,p)-ez2(m+1,n,p-1)+ez2(m,n,p-1));
            co2=coff/((zdeta(p)+zdeta(p-1))/2)/zdeta(p)*(1/kz(p)+az(p))*(1/kzm(p)+azm(p))*(ex1(m,n,p+1)-ex1(m,n,p))+coff/((zdeta(p)+zdeta(p-1))/2)/zdeta(p-1)*(1/kz(p)+az(p))*(1/kzm(p-1)+azm(p-1))*(ex1(m,n,p-1)-ex1(m,n,p))+ ex1(m,n,p);
            co3=-tdeta/e0/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(hy1(m,n,p)-hy1(m,n,p-1));
            cofff= tdeta/e0 *( pxy1(m,n,p)-bz(p)*pxz1(m,n,p));
            cofff1=-cq*(bzm(p)*myz1(m,n,p)-bxm(m)*myx1(m,n,p));
            cofff2=-cq*(bzm(p-1)*myz1(m,n,p-1)-bxm(m)*myx1(m,n,p-1));
            cofff0=-tdeta/e0/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(cofff1-cofff2);
            co=-tdeta/e0/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(ezi(m,n,p)-ezi(m,n,p-1));  %连接边界处理，位于自由空间
           so1(p-1)=co1+co2+co3+cofff0+cofff+co;
            end
            jie=matfree*so1;
            for p=2:1:pz4
            ex2(m,n,p)=jie(p-1);
            end

else  
      
            for p=2:1:pz4
            co1=tdeta/eeff(p,1)/((ydeta(n)+ydeta(n-1))/2)*(1/ky(n))*(hz1(m,n,p)-hz1(m,n-1,p))- tdeta^2/(2*u0*eeff(p,1)*xdeta(m))/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(1/kxm(m)+axm(m))*(ez2(m+1,n,p)-ez2(m,n,p)-ez2(m+1,n,p-1)+ez2(m,n,p-1));
            co2=coffe(p,1)/((zdeta(p)+zdeta(p-1))/2)/zdeta(p)*(1/kz(p)+az(p))*(1/kzm(p)+azm(p))*(ex1(m,n,p+1)-ex1(m,n,p))+coffe(p,1)/((zdeta(p)+zdeta(p-1))/2)/zdeta(p-1)*(1/kz(p)+az(p))*(1/kzm(p-1)+azm(p-1))*(ex1(m,n,p-1)-ex1(m,n,p))+De8(p,1)*ex1(m,n,p);
            co3=-tdeta/eeff(p,1)/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(hy1(m,n,p)-hy1(m,n,p-1));
            cofff= tdeta/eeff(p,1) *( pxy1(m,n,p)-bz(p)*pxz1(m,n,p));
            cofff1=-cq*(bzm(p)*myz1(m,n,p)-bxm(m)*myx1(m,n,p));
            cofff2=-cq*(bzm(p-1)*myz1(m,n,p-1)-bxm(m)*myx1(m,n,p-1));
            cofff0=-tdeta/eeff(p,1)/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(cofff1-cofff2);
co=-tdeta/e0/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(ezi(m,n,p)-ezi(m,n,p-1));  %连接边界处理
      sesan=De9(p,1)*Jpx1(m, n,p)+De10(p,1) * ex0(m, n,p)-De11(p,1)* Jpx0(m, n,p);  %色散媒质处理
           so1(p-1)=co1+co2+co3+cofff0+cofff+co +sesan;
           end
            jie=mat*so1;
            for p=2:1:pz4
            ex2(m,n,p)=jie(p-1);
            end
end
    end
end
for m=obx1:1:obx2
        for n=oby1:1:oby2
            for p=2:1:pz4
Jpx2(m,n,p)=kp2(p,1)*Jpx1(m,n,p)+kp3(p,1)*Jpx0(m,n,p)+ betaL1(p,1)* ex2(m,n,p)- betaL2(p,1)* ex0(m,n,p)+betap1(p,1)*ex1(m,n,p);
end
end
end


for  m=px1:1:px4
    for n=py1:1:py4
        for p=pz1:1:pz4+1
             mzy2(m,n,p)=bym(n)*mzy1(m,n,p)+aym(n)*(ex2(m,n+1,p)-ex2(m,n,p))/ydeta(n);
        end
    end
end

for m=px1:1:px4
    for n=py1:1:py4+1
        for p=pz1:1:pz4
            myz2(m,n,p)=bzm(p)*myz1(m,n,p)+azm(p)/2/zdeta(p)*(ex2(m,n,p+1)-ex2(m,n,p)+ex1(m,n,p+1)-ex1(m,n,p));
        end
    end
end


for m=1:1:px4
    for n=1:1:py4+1
        for p=1:1:pz4
            cofff=-cq*(bzm(p)*myz1(m,n,p)-bxm(m)*myx1(m,n,p));
            hy2(m,n,p)=cp*hy1(m,n,p)-cq*((ex1(m,n,p+1)-ex1(m,n,p)+ex2(m,n,p+1)-ex2(m,n,p))/2*(1/kzm(p)+azm(p))/zdeta(p)-(ez2(m+1,n,p)-ez2(m,n,p))*(1/kxm(m)+axm(m))/xdeta(m))+cofff+ezi(m,n,p);
         end
     end
end

for m=px1:1:px4
    for n=py1:1:py4+1
        for p=pz1+1:1:pz4
            pxz2(m,n,p)=bz(p)*pxz1(m,n,p)+az(p)*(hy2(m,n,p)-hy2(m,n,p-1)+hy1(m,n,p)-hy1(m,n,p-1))/2/((zdeta(p)+zdeta(p-1))/2);
        end
    end
end
for m=px1+1:1:px4
    for n=py1:1:py4+1
        for p=pz1:1:pz4
            pzx2(m,n,p)=bx(m)*pzx1(m,n,p)+ax(m)*(hy2(m,n,p)-hy2(m-1,n,p))/ ((xdeta(m)+xdeta(m-1))/2);
        end
    end
end

for m=i1:1:i2+1
    for p=k1:1:k2
       ezi1(m,j1-1,p)=cq*ezin2(j1-1+1)/ydeta(j1-1);
        ezi1(m,j2+1,p)=-cq/ydeta(j2+1)*ezin2(j2+1);
    end
end
 for m=i1:1:i2+1
     for n=j1:1:j2
        hxi(m,n,k1)=-tdeta/e0/2*(hxin2(n)+hxin1(n))/((zdeta(k1)+zdeta(k1-1))/2);
         hxi(m,n,k2+1)= tdeta/e0 /2*(hxin2(n)+hxin1(n))/((zdeta(k2+1)+zdeta(k2))/2);
     end
 end
for n=1:1:py4
  for m=2:1:px4
if m<obx1|m>obx2|n<oby1|n>oby2

for p=2:1:pz4
          cof2=coff/((zdeta(p)+zdeta(p-1))/2)/zdeta(p)*(1/kz(p)+az(p))*(1/kzm(p)+azm(p))*(ey1(m,n,p+1)-ey1(m,n,p))+coff/((zdeta(p)+zdeta(p-1))/2)/zdeta(p-1)*(1/kz(p)+az(p))*(1/kzm(p-1)+azm(p-1))*(-ey1(m,n,p)+ey1(m,n,p-1))+ey1(m,n,p);
          cof1=tdeta/e0/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(hx1(m,n,p)-hx1(m,n,p-1))-tdeta^2/(2*u0*e0*ydeta(n))/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(1/kym(n)+aym(n))*(ez2(m,n+1,p)-ez2(m,n,p)-ez2(m,n+1,p-1)+ez2(m,n,p-1));
          cof3=-tdeta/e0/((xdeta(m)+xdeta(m-1))/2)*(1/kx(m))*(hz1(m,n,p)-hz1(m-1,n,p));
          cofff1=-cq*(bym(n)*mxy1(m,n,p)-bzm(p)*mxz1(m,n,p));
          cofff2=-cq*(bym(n)*mxy1(m,n,p-1)-bzm(p-1)*mxz1(m,n,p-1));
          cofff0=tdeta/e0/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(cofff1-cofff2);
          cofff= tdeta/e0 *(pyz1(m,n,p)*bz(p)-pyx1(m,n,p));
          co=tdeta/e0/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(ezi1(m,n,p)-ezi1(m,n,p-1));

          so2(p-1)=cof1+cof2+cof3+cofff+cofff0+co+hxi(m,n,p);
      end
      jie=matfree*so2;
      for p=2:1:pz4
        ey2(m,n,p)=jie(p-1);
      end

else  
  
      for p=2:1:pz4
          cof2=coffe(p,1)/((zdeta(p)+zdeta(p-1))/2)/zdeta(p)*(1/kz(p)+az(p))*(1/kzm(p)+azm(p))*(ey1(m,n,p+1)-ey1(m,n,p))+coffe(p,1)/((zdeta(p)+zdeta(p-1))/2)/zdeta(p-1)*(1/kz(p)+az(p))*(1/kzm(p-1)+azm(p-1))*(-ey1(m,n,p)+ey1(m,n,p-1))+De8(p,1)*ey1(m,n,p);
          cof1=tdeta/eeff(p,1)/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(hx1(m,n,p)-hx1(m,n,p-1))-tdeta^2/(2*u0*eeff(p,1)*ydeta(n))/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))*(1/kym(n)+aym(n))*(ez2(m,n+1,p)-ez2(m,n,p)-ez2(m,n+1,p-1)+ez2(m,n,p-1));
          cof3=-tdeta/eeff(p,1)/((xdeta(m)+xdeta(m-1))/2)*(1/kx(m)+ax(m))*(hz1(m,n,p)-hz1(m-1,n,p));
          cofff1=-cq*(bym(n)*mxy1(m,n,p)-bzm(p)*mxz1(m,n,p));
          cofff2=-cq*(bym(n)*mxy1(m,n,p-1)-bzm(p-1)*mxz1(m,n,p-1));
          cofff0=tdeta/eeff(p,1)/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(cofff1-cofff2);
          cofff= tdeta/eeff(p,1) *(pyz1(m,n,p)*bz(p)-pyx1(m,n,p)*bx(m));
          co=tdeta/ eeff(p,1)/((zdeta(p)+zdeta(p-1))/2)*(1/kz(p)+az(p))/2*(ezi1(m,n,p)-ezi1(m,n,p-1));
sesan=De9(p,1)*Jpy1(m, n,p)+ De10(p,1)* ey0(m, n,p)-De11(p,1)* Jpy0(m, n,p);  %色散电流

          so2(p-1)=cof1+cof2+cof3+cofff+cofff0+co+hxi(m,n,p)+sesan;
      end
      jie=mat*so2;
      for p=2:1:pz4
        ey2(m,n,p)=jie(p-1);
      end
end
   end
end 
%
for n=oby1:1:oby2
  for m=obx1:1:obx2
      for p=2:1:pz4
Jpy2(m,n,p)=kp2(p,1)*Jpy1(m,n,p)+kp3(p,1)*Jpy0(m,n,p)+ betaL1(p,1)* ey2(m,n,p)- betaL2(p,1)* ey0(m,n,p)+betap1(p,1)*ey1(m,n,p);
end
end
end

for  m=px1:1:px4
    for n=py1:1:py4
        for p=pz1:1:pz4+1
           mzx2(m,n,p)=bxm(m)*mzx1(m,n,p)+axm(m)*(ey2(m+1,n,p)-ey2(m,n,p))/xdeta(m);
        end
    end
end
for m=px1:1:px4+1
    for n=py1:1:py4
        for p=pz1:1:pz4
         mxz2(m,n,p)=bzm(p)*mxz1(m,n,p)+azm(p)/2/zdeta(p)*(ey2(m,n,p+1)-ey2(m,n,p)+ey1(m,n,p+1)-ey1(m,n,p));
        end
    end
end

 
for m=1:1:px4+1
    for n=1:1:py4
        for p=1:1:pz4
            cofff=-cq*(bym(n)*mxy1(m,n,p)-bzm(p)*mxz1(m,n,p));
            hx2(m,n,p)=cp*hx1(m,n,p)-cq*((ez2(m,n+1,p)-ez2(m,n,p))/ydeta(n)*(1/kym(n)+aym(n))-(ey2(m,n,p+1)-ey2(m,n,p)+ey1(m,n,p+1)-ey1(m,n,p))/2/zdeta(p)*(1/kzm(p)+azm(p)))+cofff+ezi1(m,n,p);
         end
     end
end

    for m=px1:1:px4+1
    for n=py1:1:py4
        for p=pz1+1:1:pz4
            pyz2(m,n,p)=bz(p)*pyz1(m,n,p)+az(p)*(hx2(m,n,p)-hx2(m,n,p-1)+hx1(m,n,p)-hx1(m,n,p-1))/2/((zdeta(p)+zdeta(p-1))/2);
        end
    end
    end
for m=px1:1:px4+1
    for n=py1+1:1:py4
        for p=pz1:1:pz4
            pzy2(m,n,p)=by(n)*pzy1(m,n,p)+ay(n)*(hx2(m,n,p)-hx2(m,n-1,p))/ ((ydeta(n)+ydeta(n-1))/2);
        end
    end
end
for m=1:1:px4
    for n=1:1:py4
        for p=1:1:pz4+1
            hz2(m,n,p)=cp*hz1(m,n,p)-cq*((ey2(m+1,n,p)-ey2(m,n,p))/kxm(m)/xdeta(m)-(ex2(m,n+1,p)-ex2(m,n,p))/kym(n)/ydeta(n))-cq*(mzx2(m,n,p)-mzy2(m,n,p));
         end
     end
end

 for  m=px1:1:px4
    for n=py1+1:1:py4
        for p=pz1:1:pz4+1
             pxy2(m,n,p)=by(n)*pxy1(m,n,p)+ay(n)/ ((ydeta(n)+ydeta(n-1))/2)*(hz2(m,n,p)-hz2(m,n-1,p));
        end
    end
end

for  m=px1+1:1:px4
    for n=py1:1:py4
        for p=pz1:1:pz4+1
          pyx2(m,n,p)=bx(m)*pyx1(m,n,p)+ax(m)/((xdeta(m)+xdeta(m-1))/2)*(hz2(m,n,p)-hz2(m-1,n,p));
        end
    end
end


 
ex0=ex1;
ey0=ey1;
ez0=ez1;
Jpx0=Jpx1;
Jpy0=Jpy1;
Jpz0=Jpz1;
   
ex1=ex2;
ey1=ey2;
ez1=ez2;
hx1=hx2;
hy1=hy2;
hz1=hz2;


myx1=myx2;
mxy1=mxy2;
mxz1=mxz2;
mzx1=mzx2;
myz1=myz2;
mzy1=mzy2;
pyx1=pyx2;
pxy1=pxy2;
pxz1=pxz2;
pzx1=pzx2;
pyz1=pyz2;
pzy1=pzy2;
  ezin1=ezin2;
  hxin1=hxin2;
  
resulthey(t)=ey1(fieldx,fieldy,fieldz);
resulthex(t)=ex1(fieldx,fieldy,fieldz);
resulthez(t)=ez1(fieldx,fieldy,fieldz);
resulthhy(t)=hy1(fieldx,fieldy,fieldz);
resulthhx(t)=hx1(fieldx,fieldy,fieldz);
resulthhz(t)=hz1(fieldx,fieldy,fieldz);
resulthey1(t)=ey1(fieldx1,fieldy1,fieldz1);
resulthex1(t)=ex1(fieldx1,fieldy1,fieldz1);
resulthez1(t)=ez1(fieldx1,fieldy1,fieldz1);
resulthhy1(t)=hy1(fieldx1,fieldy1,fieldz1);
resulthhx1(t)=hx1(fieldx1,fieldy1,fieldz1);
resulthhz1(t)=hz1(fieldx1,fieldy1,fieldz1);
resulthey2(t)=ey1(fieldx2,fieldy2,fieldz2);
resulthex2(t)=ex1(fieldx2,fieldy2,fieldz2);
resulthez2(t)=ez1(fieldx2,fieldy2,fieldz2);
resulthhy2(t)=hy1(fieldx2,fieldy2,fieldz2);
resulthhx2(t)=hx1(fieldx2,fieldy2,fieldz2);
resulthhz2(t)=hz1(fieldx2,fieldy2,fieldz2);

end
%plot(xzh,resulthez2);
plot(xzh,resulthez1);
%%%%%
toc
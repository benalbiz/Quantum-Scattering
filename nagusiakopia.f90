program nagusia
use tipoak
use gauss_legendre

integer, parameter:: ck=selected_real_kind(16)
integer,parameter:: elxi=30,elphi=30
integer, parameter::el=elxi*elphi,ezmug=4*elphi,ezbar=16*elphi*elxi
integer,parameter:: hdim=3*(8*elphi+12)+1,ezt=ezmug+ezbar, digkanp=8*elphi+12,ldim=3*digkanp+1
integer,parameter::ngaus=8


integer:: i,ipp,kr,kp


!!!elr: r-ko elementu kant. elphi: phi-ko elementu kant; el: el. kant totala.
!ezmug: mugaldeko ezezagun kop. ezbar: barneko ezezagun kop. hdim: matrize globalaren dimentsioak

complex(kind=ck), parameter:: irud=(0.0,1.0)
real(kind=dp), parameter:: pi=acos(-1.0_dp)
real(kind=dp):: phi,phi0,phil,xi0,xil,xi,ximax,k,a,amult,rmult,sekzint,sekz
real(kind=dp)::hxi,hphi,rg,pg,rg2,vol,u,lerro,zutabe,der,angelua
real(kind=dp),dimension(36,36):: hl
real(kind=dp),dimension(6)::wr,wp
real(kind=dp),dimension(6,ngaus)::p,pd,pdd
real(kind=dp),dimension(0:2*elphi):: sigma
!hl matrize lokala
      !xi=k*r; alpha=a*k, non a zilindroaren erradioa den; sigma sekzio
      !eragilea; phi angelua; r koordenatu erradiala.


complex(kind=ck),dimension(ldim,ezt):: h !matrize globala
complex(kind=ck),dimension(ezt):: bek, ipiv,zerobek !emaitzen bektorea eta bektore huts
real(kind=dp),dimension(ezt):: bekimag,bekre,bekam


real(kind=dp),dimension(el,36):: mlog
complex(kind=ck)::usar,uirt,usarxi,usarphi,uirtxi,usarxiphi,uirtxiphi
!errad. zonaldean, psi=uhinlau+f(phi)*uhinerro
!uhinlau=exp(i*cos(phi)*r)
!uhinerro=exp(i*r)*exp(pi*i/4)/sqrt(r) 

real(kind=dp),dimension(ngaus)::root,alpha
integer::ixi,iphi,n,ixin,jpun,ixil,jphil,iphin,ixipun,iphipun,j
integer:: jphin,jxin,jxipun,jphipun,ipun,iphil,il,jl,igb,ig,jg,jxil,info

mlog=matglob(elxi,elphi,el)


!!!barreraren erradioa eta erradio maximoa:
ximax=pi*10.0_dp
a=pi
!!!
hxi=ximax/(2*elxi)
hphi=pi/(2*elphi)


open(unit=15,file="uhinfuntzioa.dat",action="write",status="replace")
open(unit=16,file="sekzioa.dat",action="write",status="replace")
open(unit=17,file="xi0.dat",action="write",status="replace")
bek=0.0
zerobek=0.0
hl=0.0
h=0.0 

!polinomioen integral zuzen



 do i=1,6,2
     wr(i)=1.0_dp
     wp(i)=1.0_dp
 end do
   do i=2,6,2
     wr(i)=hxi
     wp(i)=hphi
  end do

 call gauss(ngaus,alpha,root)

 call poly(ngaus,root,p,pd,pdd)

!ekuaziorako: Pj/h*d^2Pi/d^2xil+Pj*(xi0+hxil)/h*dPi/dxil+Pi/h*dPj/dphil+(1-u)PiPj=0

!indize globaletatik hasita, elementuz-elementu egin behar da kalkulua. ondo ordenatu behar dira!!
!indize globaletan eginda, jada jarraitua izango da funtzioa.


!do elementuetan
do ixi=1,elxi
do iphi=1,elphi
      n=elxi*(ixi-1)+iphi !elementuaren indizea

!do elementu barruan


!do i indizeak

do ixin=1,3 !nodo bakoitzaren xi indize
do iphin=1,3 !nodo bakoitzaren phi indizea
     do ixipun=0,1 !nodo baten barruko koefiziente baten xi indizea
     do  iphipun=0,1 !nodo baten barruko koefiziente baten phi indizea
     ipun=ixipun*2+iphipun+1 !nodo baten barruko koefiziente bakoitzaren indize totala
     ixil=(ixin-1)*2+ixipun+1 !nodo bakoitzerako polinomioen oinarrian xi indizea (1-6)
     iphil=(iphin-1)*2+iphipun+1 !nodo bakoitzerako polinomioen oinarrian phi indizea (1-6)

     il=(ixin-1)*12+(iphin-1)*4+ipun !koefiziente bakoitzaren indize lokala
     ig=mlog(n,il)
    
     if (ig==0) then
             cycle
     end if

!do j indizeak

do jxin=1,3
do jphin=1,3
     do jxipun=0,1
     do jphipun=0,1
     jpun=jxipun*2+jphipun+1
     jxil=(jxin-1)*2+jxipun+1
     jphil=(jphin-1)*2+jphipun+1

     jl=(jxin-1)*12+(jphin-1)*4+jpun
     jg=mlog(n,il)

     if (jg==0) then
             cycle
     end if

!zero matrize lokala

do  kr=1,ngaus
              do  kp=1,ngaus
              rg=(ixi-1)*2*hxi+(root(kr)+1)*hxi
              pg=(iphi-1)*2*hphi+(root(kp)+1)*hphi
              rg2=rg*rg
              vol=alpha(kr)*alpha(kp)*rg*hxi*hphi*wr(ixil)*wp(iphil)*wr(jxil)*wp(jphil)
              

u=0.0_dp
if (xi<=a) then
        u=1000000.0_dp
end if

lerro= p(ixil,kr) *p(iphil,kp)*vol
zutabe= p(jxil,kr)*p(jphil,kp)
der= pdd(jxil,kr)/(hxi**2)*p(jphil,kp)+pd(jxil,kr)/rg/hxi*p(jphil,kp)+p(jxil,kr)*pdd(jphil,kp)/(hphi**2)/rg2 



hl(il,jl)=hl(il,jl)+lerro*(der+zutabe*(1.0_dp-u))


end do
end do




!indize globala/banded storage????
igb=ig-jg+(2*digkanp+1)


!!! xi=ximax mugalde baldintzak

if (ixi==elxi .AND. jxin==3) then
        angelua=(iphi-1)*hphi*2+(jphin-1)*hphi
        
end if

usar=exp(irud*ximax*cos(angelua))
uirt=exp(irud*ximax)/sqrt(ximax)*exp(irud*pi/4)

if (jpun==1) then
        h(igb,jg)=h(igb,jg)+hl(il,il)*uirt
        bek(ig)=bek(ig)-hl(il,jl)*usar
else if (jpun==2) then
        usarphi=-irud*ximax*sin(angelua)*usar !dusar/dphi
        h(igb,jg)=h(igb,jg)+hl(il,jl)*uirt
        bek(ig)=bek(ig)-hl(il,jl)*usarphi
else if (jpun==3) then
        usarxi=irud*cos(angelua)*usar !dusar/dxi
        uirtxi=(irud-0.5_dp/rg)*uirt !duirt/dxi
        h(igb,jg)=h(igb,jg)+hl(il,jl)*uirtxi
        bek(ig)=bek(ifg)-hl(il,jl)*usarxi
else if (jpun==4) then
        usarxiphi=(-irud*sin(angelua)+ximax*cos(angelua)*sin(angelua))*usar
        uirtxiphi=(irud-0.5_dp/rg)*uirt
        h(igb,jg)=h(igb,jg)+hl(il,jl)*uirtxiphi
        bek(ig)=bek(ig)-hl(il,jl)*usarxiphi
else
        h(igb,jg)=h(igb,jg)+hl(il,jl)
end if

end do 
end do 
end do 
end do 
end do
end do 
end do 
end do 
end do 
end do 


call zgbsv(ezt,digkanp,digkanp,1,h,ldim,ipiv,bek,ezt,info)

bekre=real(bek)

bekimag=aimag(bek)



write(unit=15,fmt="(2f15.4)") bek




!!! sekzio eragilea teorema optikoarekin

sekz=2*sqrt(2*pi)*aimag(bek(ezbar+1))



!bektorearen azkenengo gaiak, sigma, dsigma/dphi

sigma(0)=abs(bek(ezbar+1))**2

do i=1,2*elphi
sigma(i)=abs(bek(ezbar+i*2))**2
end do

!sekzint: sekzio eragilea integratzen

sekzint=0.0_dp
do i=-2*elphi,2*elphi
ipp=i
if (i<=0) then
        ipp=-i
end if
write(unit=16,fmt="(2f12.6)") i*hphi, sigma(ipp)

if (i<=0) then
        cycle
end if
sekzint=sekzint+(sigma(i)+sigma(i-1))*pi/(2*elphi)

end do 

close(unit=15)
close(unit=16)
close(unit=17)
!!!!!!! Erabiltzeko LAPACK honela konpilatu: F zuenprogrma.f90 -llapack -lblas

!!!!!SUBRUTINA: call zgsbv(ETC)... 
contains

    
     function matglob(elxi,elphi,el)
     integer,intent(in):: elxi, elphi, el
     real(kind=dp),dimension(el,36)::mc,matglob
     
     integer:: nlerro

     integer:: i,j,k,l,m,n,ixi,iphi,nel,ipun,il,ezbar

      nlerro=elphi*8
       mc=0.0_dp


     !elementuen gaineko do
     do i=1,elxi  !i: r-ko elementu zenbakia
     do j=1,3    !j: r-ko indize lokala
     ixi=(i-1)*2+j   !ir: phi-ko indize globala
     do k=1,elphi  !k: phi-ko elementu zenbakia
     do l=1,3     !l: phi-ko indize lokala
     iphi=(k-1)*2+l !!! phi-ko indize globala

     nel=(i-1)*elphi+k !!!elementuaren indizea

     !!nodoen gaineko do
     do m=0,1
     do n=0,1
     ipun=m*2+n+1 !nodo batean, koefiziente baten indizea
     il=(j-1)*12+(l-1)*4+ipun !koefiziente baten indize lokala

     !! mugalde baldintzak xi=ximax-en

     if (ixi==2*elxi+1) then
           ezbar=nlerro*elxi*2
           if (iphi==1) then
             if (n==1) then
                    mc(nel,il)=0
             end if
             if (n==0) then
                    mc(nel,il)=ezbar+1
             end if
          else if (iphi==2*elphi+1) then
             if (n==1) then
                     mc(nel,il)=0
             end if
             if (n==0) then
                    mc(nel,il)=ezbar+1+(iphi-2)*2+1
             end if
          else
             mc(nel,il)=ezbar+1+(iphi-2)*2+1+n
          end if
          exit

     end if

     !! dpsi/dphi=0 phi=0, phi=pi

     if (iphi==1) then
             if (n==1) then
                    mc(nel,il)=0
             end if
             if (n==0) then
                    mc(nel,il)=(ixi-1)*nlerro+m+1
             end if
     else if (iphi==2*elphi+1) then
             if (n==1) then
                     mc(nel,il)=0
             end if
             if (n==0) then
                    mc(nel,il)=(ixi-1)*nlerro+2+(iphi-2)*4+m+1
             end if
     else
             mc(nel,il)=(ixi-1)*nlerro+2+(iphi-2)*4+ipun
     end if

     cycle

     end do
     end do
     end do
     end do
     end do
     end do 
    

     matglob=mc
     end function matglob     


 subroutine poly(ngaus,root,phi,phid,phidd)
    integer,intent(in)::ngaus
    real(kind=dp),dimension(:),intent(in):: root


    real(kind=dp),dimension(:,:),intent(inout):: phi,phid,phidd
    real(kind=dp)::x,x0,x1,x2,x3,x4,x5
    integer::n    
         do n = 1, ngaus
                x=root(n)
                x0=1.0_dp
                x1=x
                x2=x*x
                x3=x2*x
                x4=x3*x
                x5=x4*x

 phi  (1,n) =              x2   -  5*x3/4 -    x4/2 +  3*x5/4
 phid (1,n) =            2*x1   - 15*x2/4 -  2*x3   + 15*x4/4
 phidd(1,n) =            2*x0   - 15*x1/2 -  6*x2   + 15*x3

 phi  (2,n) =              x2/4 -    x3/4 -    x4/4 +    x5/4
 phid (2,n) =              x1/2 -  3*x2/4 -    x3   +  5*x4/4
 phidd(2,n) =              x0/2 -  3*x1/2 -  3*x2   +  5*x3

 phi  (3,n) = x0       - 2*x2             +    x4
 phid (3,n) =          - 4*x1             +  4*x3
 phidd(3,n) =          - 4*x0             + 12*x2

 phi  (4,n) =      x1           -  2*x3             +    x5
 phid (4,n) =      x0           -  6*x2             +  5*x4
 phidd(4,n) =                   - 12*x1             + 20*x3

 phi  (5,n) =              x2   +  5*x3/4 -    x4/2 -  3*x5/4
 phid (5,n) =            2*x1   + 15*x2/4 -  2*x3   - 15*x4/4
 phidd(5,n) =            2*x0   + 15*x1/2 -  6*x2   - 15*x3

 phi  (6,n) =            - x2/4 -    x3/4 +    x4/4 +    x5/4
 phid (6,n) =            - x1/2 -  3*x2/4 +    x3   +  5*x4/4 
 phidd(6,n) =            - x0/2 -  3*x1/2 +  3*x2   +  5*x3

  end do
  end subroutine poly

      subroutine gauss (NGAUS,ALPHA,ROOT)
      integer,intent(in):: ngaus
      real(kind=dp),dimension(:),intent(inout):: ALPHA, ROOT
      integer::imax,i
! Modfiy for higher quadrature if needed

      
      ALPHA(1)=0.362683783378362_dp
      ALPHA(2)=0.313706645877887_dp
      ALPHA(3)=0.222381034453374_dp
      ALPHA(4)=0.101228536290376_dp
      ROOT(1)=0.183434642495650_dp
      ROOT(2)=0.525532409916329_dp
      ROOT(3)=0.796666477413627_dp
      ROOT(4)=0.960289856497536_dp
 
      imax=ngaus/2
      DO  i=1,imax
      ALPHA(NGAUS+1-I)=ALPHA(I)
      ROOT(NGAUS+1-I)=-ROOT(I)
      end do
       
      end subroutine gauss



end program nagusia

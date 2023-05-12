program nagusia
use tipoak
use polinomioak
use gauss_legendre

integer, parameter:: ck=selected_real_kind(16)
integer,parameter:: elxi=30,elphi=30
integer, parameter::el=elxi*elphi,ezmug=4*elphi,ezbar=16*elphi*elxi
integer,parameter:: hdim=3*(8*elphi+12)+1,ezt=ezmug+ezbar, digkanp=8*elphi+12,ldim=3*digkanp+1

integer:: i,ipp


!!!elr: r-ko elementu kant. elphi: phi-ko elementu kant; el: el. kant totala.
!ezmug: mugaldeko ezezagun kop. ezbar: barneko ezezagun kop. hdim: matrize globalaren dimentsioak

complex(kind=ck), parameter:: irud=(0.0,1.0)
real(kind=dp), parameter:: pi=acos(-1.0_dp)
real(kind=dp):: phi,phi0,phil,xi0,xil,alpha,xi,ximax,k,a,amult,rmult,sekzint,sekz
real(kind=dp)::hxi,hphi
real(kind=dp),dimension(36,36):: hl
real(kind=dp),dimension(6)::pint,bint,cint,derbint
real(kind=dp),dimension(0:2*elphi):: sigma
!hl matrize lokala
      !xi=k*r; alpha=a*k, non a zilindroaren erradioa den; sigma sekzio
      !eragilea; phi angelua; r koordenatu erradiala.


complex(kind=ck),dimension(ldim,ezt):: h !matrize globala
complex(kind=ck),dimension(ezt):: bek, ipiv,zerobek !emaitzen bektorea eta bektore huts

real(kind=dp),dimension(el,36):: mlog
complex(kind=ck)::usar,uirt,usarxi,usarphi,uirtxi,usarxiphi,uirtxiphi
!errad. zonaldean, psi=uhinlau+f(phi)*uhinerro
!uhinlau=exp(i*cos(phi)*r)
!uhinerro=exp(i*r)*exp(pi*i/4)/sqrt(r) 


mlog=matglob(elr,elphi,el)


!!!barreraren erradioa eta erradio maximoa:
ximax=pi*20
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

!polinomioen integral zuzena

pint(1)=integ(p1)
pint(2)=integ(p2)
pint(3)=integ(p3)
pint(4)=integ(p4)
pint(5)=integ(p5)
pint(6)=integ(p6)

derbint(1)=integ(derb1)
derbint(2)=integ(derb2)
derbint(3)=integ(derb3)
derbint(4)=integ(derb4)
derbint(5)=integ(derb5)
derbint(6)=integ(derb6)



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
     ig=matglob(n,il)
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
     jg=matglob(n,il)

     if (jg==0) then
             cycle
     end if

!zero matrize lokala


xi0=ixi*hxi
xi=(ixi-1)*hxi*2+hxi*(ixin-1)
! integrazioa (i,j): pint(j)*derbint(i)/hxi^2+pint(j)*bint(i)+pint(i)*cint(j)+(1-u)*pint(i)*pint(j)

write(unit=17, fmt="(3f12.12)") xi0,hxi,hphi


u=0.0_dp
if (xi<=a) then
        u=1000000.0_dp
end if

!j-rekiko integrazioa:
batugai=(1-u)*pint(jxil)*pint(jphil)


if (jxil==1) then
        batugai=batugai+pint(jphil)*derbint(jxil)/hxi**2+pint(jphil)*integ(b1)+pint(jxil)*integ(c1)
else if (jxil==2) then
        batugai=batugai+pint(jphil)*derbint(jxil)/hxi**2+pint(jphil)*integ(b2)+pint(jxil)*integ(c2)
else if (jxil==3) then
        batugai=batugai+pint(jphil)*derbint(jxil)/hxi**2+pint(jphil)*integ(b3)+pint(jxil)*integ(c3)
else if (jxil==4) then
        batugai=batugai+pint(jphil)*derbint(jxil)/hxi**2+pint(jphil)*integ(b4)+pint(jxil)*integ(c4)
else if (jxil==5) then
        batugai=batugai+pint(jphil)*derbint(jxil)/hxi**2+pint(jphil)*integ(b5)+pint(jxil)*integ(c5)
else 
        batugai=batugai+pint(jphil)*derbint(jxil)/hxi**2+pint(jphil)*integ(b6)+pint(jxil)*integ(c6)
end if

!ez dakit zelan sartu i indizearen menpekotasuna

hl(il,jl)=hl(il,jl)+batugai




!!!!MUGALDE BALDINTZA PERIODIKOAK: -1=N
!!!!ZUZENEAN PSI(XI<ALFA)=0


!!probatzeko
if (il==1) then
        print*, hl

end if


!indize globala/banded storage????
igb=ig-jg+(2*digkanp+1)


!!! xi=ximax mugalde baldintzak

if (ixi==elxi .AND. jxin==3) then
        angelua=(elphi-1)*hphi*2+(jphin-1)*hphi
end if

usar=zexp(i*ximax*cos(phi))
uirt=zexp(i*ximax)/sqrt(ximax)*zexp(i*pi/4)

if (jpun==1) then
        h(igb,jg)=h(igb,jg)+hl(il,il)*uhinirt
        bek(ig)=bek(ig)-hl(il,jl)*uhinsar
else if (jpun==2) then
        usarphi=-irud*ximax*sin(angelua)*usar
        h(igb,jg)=h(igb,jg)+hl(il,jl)*uirt
        bek(ig)=bek(ig)-hl(il,jl)*usarphi
else if (jpun==3) then
        usarxi=irud*cos(angelua)*usar
        uirtxi=(irud-0.5_dp/xi)*uirt
        h(igb,jg)=h(igb,jg)+hl(il,jl)*uirtxi
        bek(ig)=bek(ig)-hl(il,jl)*usarxi
else if (jpun==4) then
        usarxiphi=(-irud*sin(angelua)+ximax*cos(angelua)*sin(angelua))*usar
        uirtxiphi=(irud-0.5_dp/xi)*uirt
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


call zgsbv(ezt,digkanp,digkanp,1,h,ldim,ipiv,bek,ezt,info)

write(unit=15,fmt="(2d12.4)") bek

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
sekzint=sekzint+(f(i)+f(i-1))*pi/(2*elphi)

end do 

close(unit=15)
close(unit=16)
close(unit=17)
!!!!!!! Erabiltzeko LAPACK honela konpilatu: F zuenprogrma.f90 -llapack -lblas

!!!!!SUBRUTINA: call zgsbv(ETC)... 
contains

    
     function matglob(elxi,elphi,el)
     integer,intent(in):: elxi, elphi, el
     real(kind=dp),dimension(el,36):: matglob,mc
     
     integer:: nlerro

     integer:: i,j,k,l,m,n,ixi,iphi,nel,ipun,il,ezbar

      nlerro=elphi*8
       mc=0.0


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
             mc(nel,il)=ezbar+1+(iphi-2)*2+1
          end if
          cycle

     end if

     !! dpsi/dphi=0 phi=0, phi=pi

     if (iphi==1) then
             if (n==1) then
                    mc(nel,il)=0
             end if
             if (n==0) then
                    mc(nel,il)=(j-1)*nlerro+m+1
             end if
     else if (iphi==2*elphi+1) then
             if (n==1) then
                     mc(nel,il)=0
             end if
             if (n==0) then
                    mc(nel,il)=(j-1)*nlerro+2+(iphi-2)*4+m+1
             end if
     else
             mc(nel,il)=(j-1)*nlerro+2+(iphi-2)*4+ipun
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


    function b1(x)
            real(kind=dp),intent(in)::x
            real(kind=dp)::b1,xi0,hxi
            read(unit=17,fmt=*) xi0,hxi
            b1=bgaia1(xi0,hxi,x)
       end function b1

       function b2(x)
            real(kind=dp),intent(in)::x
            real(kind=dp)::b2,xi0,hxi
            read(unit=17,fmt=*) xi0,hxi

            b2=bgaia2(xi0,hxi,x)
       end function b2

       function b3(x)
            real(kind=dp),intent(in)::x
            real(kind=dp)::b3,xi0,hxi
            read(unit=17,fmt=*) xi0,hxi
            b3=bgaia3(xi0,hxi,x)
       end function b3

       function b4(x)
             real(kind=dp),intent(in)::x
             real(kind=dp)::b4,xi0,hxi
             read(unit=17,fmt=*) xi0,hxi

             b4=bgaia4(xi0,hxi,x)
        end function b4

        function b5(x)
             real(kind=dp),intent(in)::x
             real(kind=dp)::b5
             b5=bgaia5(xi0,hxi,x)
         end function b5

         function b6(x)
             real(kind=dp),intent(in)::x
             real(kind=dp)::b6,xi0,hxi
             read(unit=17,fmt=*) xi0,hxi
             b6=bgaia6(xi0,hxi,x)
         end function b6
          
         function c1(x)
             real(kind=dp),intent(in)::x
             real(kind=dp)::c1,xi0,hxi,hphi
             read(unit=17,fmt=*) xi0,hxi,hphi
             c1=p1dd(xi0,hxi,hphi,x)
         end function c1
         
         function c2(x)
              real(kind=dp),intent(in)::x
              real(kind=dp)::c2,xi0,hxi,hphi
              read(unit=17,fmt=*) xi0,hxi,hphi
              c2=p2dd(xi0,hxi,hphi,x)
          end function c2

          function c3(x)
               real(kind=dp),intent(in)::x
               real(kind=dp)::c3,xi0,hxi,hphi
               read(unit=17,fmt=*) xi0,hxi,hphi
               c3=p3dd(xi0,hxi,hphi,x)
          end function c3

          function c4(x)
               real(kind=dp),intent(in)::x
               real(kind=dp)::c4,xi0,hxi,hphi
               read(unit=17,fmt=*) xi0,hxi,hphi
               c4=p4dd(xi0,hxi,hphi,x)
           end function c4
          
           function c5(x)
               real(kind=dp),intent(in)::x
               real(kind=dp)::c5,xi0,hxi,hphi
               read(unit=17,fmt=*) xi0,hxi,hphi
               c5=p5dd(xi0,hxi,hphi,x)
           end function c5

           function c6(x)
               real(kind=dp),intent(in)::x
               real(kind=dp)::c6,xi0,hxi,hphi
               read(unit=17,fmt=*) xi0,hxi,hphi
               c6=p6dd(xi0,hxi,hphi,x)
           end function c6

                function integ(g,x)
        real(kind=dp),intent(in)::x
        real(kind=dp)::a,b,integ,em
        interface

        function g(x)
               real,intent(in)::x
               real::g
        end function g
        end interface

        a=-1.0_dp
        b=1.0_dp
        call qgauss(g,a,b,em,10)
        integ=em
        end function integ



end program nagusia

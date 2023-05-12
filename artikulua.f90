! Janine Shertzer   October 1, 2015 

!     FEM code to solve the problem of particles  scattering  off of 
!     an infinte cylindrical barrier (circular, square or diamond) 
!     using dimensionless cylindrical coordinates r=rho*k  p=phi 

! range p [0:pi]  require dpsi/dp=0 at p=0 and p=pi
! range r [0:rmax] require psi=inc+out*f at r=rmax

!  NELP   number of elements in p
!  NELR   number of elements in r
!  NEL    total number of elements

!  NGRID  number of unknowns on the interior grid points  
!    psi, dpsi/dp, dpsi/dr, d^2psi/drdp (except at p=0,pi, where dpsi/dp=0)
!    [2+4*(2*nelp-1)+2]*(2*nelr)=16*nelp*nelr   

!  NBC   number of unknowns on the boundary  (f, df/dp)
!    1+2*(2*nelp-1)+1=4*nelp           
!  NNN = NGRID+NBC total number of unknowns

!  NOD = 8*NELP+12  number of off-diagonal elements in the global matrix
!  LDIM = 3*NOD+1 leading dimension of global matrix (includes workspace)
!  NGAUS = 8  Gauss quadrature for numerical integration 
!  H(LDIM,NNN)  global matrix in banded storage mode

!  Lapack routine zgbsv (and dependencies) to solve Hx=y
!  VEC(NNN)    VEC=Y on entry to zgbsv,  VEC=X on exit  
!  IPIV(NNN)   working space

!  PHI(6,NGAUS), PHID(6,NGAUS), PHIDD(6,NGAUS) 6 FEM basis functions 
!  and derivatives evaluated at the quadrature points
!  In the paper, these functions are referred to as P_i(x)
!  WR(6), WP(6) scaling factors for FEM basis functions
!
!  ALPHA(NGAUS), ROOT(NGAUS)   weights and roots for Gauss quadrature
!  HL(36,36),UL(36,36)  local matrices, generated for each element
!  MC(NEL,36)   connectivity matrix which maps local index to global index
!  F(0:2*NELP) f and df/dp  on boundary

!  Code is currently set to produce the first line in the convergence table.


program scatter
  IMPLICIT REAL*8 (a-h,o-z)
  parameter (nelr= 30, nelp= 30,nel = nelr*nelp)
  parameter (ngrid=nelr*nelp*16, nbc=nelp*4)
  parameter (nnn = ngrid+nbc,nod=nelp*16+12,ldim=3*nod+1)
  parameter (ngaus = 8)
  parameter (job=1)
  dimension phi(6,ngaus), phid(6,ngaus),phidd(6,ngaus),alpha(ngaus),root(ngaus) 
  dimension hl(36,36),ul(36,36),mc(nel,36),wr(6),wp(6),f(0:2*nelp)
  dimension ipiv(nnn)
  dimension emaitza(3)
   

  complex*16 h(ldim,nnn),vec(nnn)
  complex*16 one,imag,zero,psi,psiem
  complex*16 inc,incr,incp,incpr,out,outr,outp,outpr

  data one   / (1.0D0, 0.0D0) /
  data imag  / (0.0D0, 1.0D0) /
  data zero  / (0.0D0, 0.0D0) /
  data pi/3.141592654d0/
  integer l

  


! data files for scattering amplitude and wave function
  open (unit=18, file='SEKZIOA.amp')
  open (unit=19, file='KOEFIZIENTEAK.wave')
  open (unit=20, file="EMAITZA.dat")


! barrier radius a  and  cutoff rmax 
  amult= 1.d0
  a= pi*amult
  rmult= 10.0d0
  rmax=pi*rmult


! calculate total storage
  stor=16.d-9*ldim*nnn

! print out parameters to screen and file
  write (18,1001) amult, rmult,ngaus
  write ( 6,1001) amult,rmult,ngaus
1001   format ('#  a = ',f5.2,' pi     r_max = ',f5.2, ' pi    ngaus  =',i3)
  write (18,1002) nelr,nelp,nnn,stor
  write ( 6,1002) nelr,nelp,nnn,stor    
1002  format ('#  nelr = ',i4,'     nelp = ',i4,/    '#  nnn = ',i8,'   stor = ',f6.2, 'Gbytes'/)

! set up weight factors for FEM basis functions 
  hr=rmax/nelr/2
  hr2=hr*hr
  hp=pi/nelp/2
  hp2=hp*hp

   do 10 i=1,6,2
     wr(i)=1.0d0
     wp(i)=1.0d0
10   continue
   do 20 i=2,6,2
     wr(i)=hr
     wp(i)=hp
20   continue

!  Gauss quadrature weights and roots
!  currently using ngaus=8
   call gauss(ngaus,alpha,root)

!  create mc which maps local to global index
   call conmat(nelp,nelr,nel,mc)

!  FEM basis functions and derivatives
   call poly(ngaus,root,phi,phid,phidd)
  
!  zero global matrix and vector
   h=0.0d0
   vec=0.0d0

! do loop ever elements
    do 400 ielr=1,nelr
    do 390 ielp=1,nelp
    iel=(ielr-1)*nelp+ielp

! do loop for local index il
      do 300 irn=1,3
      do 290 ipn=1,3
      do 280 ird=0,1
      do 270 ipd=0,1
      id=ird*2+ipd+1
      irl=(irn-1)*2+ird+1
      ipl=(ipn-1)*2+ipd+1

      il=(irn-1)*12+(ipn-1)*4+id
      ig=mc(iel,il)
      if (ig.eq.0) go to 270

! do loop for local index jl
         do 200 jrn=1,3
         do 190 jpn=1,3
         do 180 jrd=0,1
         do 170 jpd=0,1
         jd=jrd*2+jpd+1
         jrl=(jrn-1)*2+jrd+1
         jpl=(jpn-1)*2+jpd+1
         jl=(jrn-1)*12+(jpn-1)*4+jd
         jg=mc(iel,jl)
         if (jg.eq.0) go to 170

! zero local matrix
         ul=0.0d0
         hl=0.0d0   

! do loop for numerical integration 
              do 100 kr=1,ngaus
              do  90 kp=1,ngaus
              rg=(ielr-1)*2*hr+(root(kr)+1)*hr
              pg=(ielp-1)*2*hp+(root(kp)+1)*hp
              rg2=rg*rg
              vol=alpha(kr)*alpha(kp)*rg*hr*hp*wr(irl)*wp(ipl)*wr(jrl)*wp(jpl)
              pot=0.0d0
!  POTENTIAL
!  for diamond barrier, replace pg in next two lines with pg+pi/4
              x=rg*dcos(pg)
              y=rg*dsin(pg)

! cylindrical barrier
     !  if (rg.lt.a) pot=1.d6
!  square or diamond barrier 
   if ((dabs(x).le.a).and.(dabs(y).le.a)) pot=1.d6

              row= phi  (irl,kr)      *phi(ipl,kp)*vol
              col= phi  (jrl,kr)      *phi(jpl,kp)
              der= phidd(jrl,kr)/hr2  *phi(jpl,kp) &
                  +phid (jrl,kr)/rg/hr*phi(jpl,kp) &
                  +phi(jrl,kr)      *phidd(jpl,kp)/hp2/rg2 
!  local matrix element 
              hl(il,jl)=hl(il,jl)+row*(der+col*(1.0d0-pot))
              90   continue
             100   continue

! igb - global index in banded storage   
       igb=ig-jg+(2*nod+1)

!  scattering boundary conditions at r=rmax
       if ((ielr.eq.nelr).and.(jrn.eq.3)) then
             angle=(ielp-1)*hp*2+(jpn-1)*hp

!  inc = incident wave  exp(i*r*cos(p))
!  out = outgoing wave  exp(i*r)/sqrt(r)*exp(i*pi/4)
!  add on contribution to global matrix h and vector  vec

             inc=zexp(imag*rmax*dcos(angle))
             out=zexp(imag*rmax)/dsqrt(rmax)*zexp(imag*pi/4)
               if (jd.eq.1) then
                  h(igb,jg)=h(igb,jg)+hl(il,jl)*out
                  vec(ig) =vec(ig) -hl(il,jl)*inc
               else  if (jd.eq.2) then
                  incp=-imag*rmax*dsin(angle)*inc
                  outp=out
                  h(igb,jg)=h(igb,jg)+hl(il,jl)*outp
                  vec(ig) =vec(ig) -hl(il,jl)*incp
               else if (jd.eq.3) then
                   incr=imag*dcos(angle)*inc
                   outr=(imag-.5d0/rg)*out
                   h(igb,jg)=h(igb,jg)+hl(il,jl)*outr
                   vec(ig) =vec(ig) -hl(il,jl)*incr
                else  if (jd.eq.4) then
                   incpr=(-imag*dsin(angle)+rmax*dcos(angle)*dsin(angle))*inc
                   outpr=(imag-.5/rg)*out
                   h(igb,jg)=h(igb,jg)+hl(il,jl)*outpr
                   vec(ig) =vec(ig) -hl(il,jl)*incpr           
             endif 
       else 

!  interior grid points - add local matrix to global matrix
       h(igb,jg)=h(igb,jg)+hl(il,jl)
       endif
         
170   continue
180   continue
190   continue
200   continue
270   continue
280   continue
290   continue
300   continue
390   continue
400   continue

!  Solve linear equations using lapack routine (netlib.org)

     call zgbsv(nnn,nod,nod,job,h,ldim,ipiv,vec,nnn,info)
     
     write (19,1007) vec
1007 format (2d12.4)

!grafikatzeko

angelua=0
angelua0=0
l=0
do i=1,nel

xi0=int(i/nelp)*2*hr

angelua0=(i-1)*hp*2-int(i/nelp)*nelp*2*hp

psiem=0
do j=1,3
   xi=xi0+(j-1)*hr
      
   xil=real(j-2,8)
   
   do k=1,3
   angelua=angelua0+(k-1)*hp
    
   phil=real(k-2,8)
      call psigaia(i,xil,phil,vec,psiem)
      
      
      l=l+1
      psi2=abs(psiem)**2
      
      write(20,3000) xi*cos(angelua),xi*sin(angelua),psi2
      3000 format(30f14.8)
    end do

end do

end do

print*, l, ngrid
!gehigarria


! cross1 - total cross section using optical theorem
    cross1=2*dsqrt(2*pi)*aimag(vec(ngrid+1))

! last entries in vec are value of f and df/dphi 
    f(0)=zabs(vec(ngrid+1))**2
    do 410 i=1,2*nelp
    f(i)=zabs(vec(ngrid+i*2))**2
410 continue

!  cross2 - total cross section by integrating f
    cross2=0.0d0
    do 500 ip=-2*nelp,2*nelp
    ipp=ip
    if (ip.le.0) ipp=-ip
    write (18,18) ip*hp, f(ipp)
18  format (2f12.6)
    if (ip.le.0) go to 500
    cross2=cross2+(f(ip)+f(ip-1))*pi/(2*nelp)
500  continue

!  (OPTIONAL) Print out values of f(phi)^2 for convergence study
      write (6,1010)  f(0)
      write (6,1011)  f(2*nelp/10)
      write (6,1012)  f(nelp)
      write (6,1013) f(2*nelp)
1010   format (' |f(0)|^2     = ',f8.3)
1011   format (' |f(pi/10)|^2 = ',f8.3)
1012   format (' |f(pi/2)|^2  = ',f8.3)
1013   format (' |f(pi)|^2    = ',f8.3)

     write (6,1006) cross1, cross2
     write (18,1006) cross1, cross2
1006  format ('#  sigma_opt = ',f8.3,'  sigma_int = ',f8.3)

stop
contains

!_________________________________________________________________________________

!    CONMAT generates mc(nel,36) which maps local index il to the global index ig

          subroutine conmat (nelp,nelr,nel,mc)
          implicit real*8 (a-h,o-z)
          dimension mc(nel,36)
          nrow=nelp*8
          igrid=nrow*(nelr*2) !!!hau ez dabil ez dakit zeba
          mc=0
! do loop over elements
          do 100 ielr=1,nelr
          do 90 irn=1,3
          ir=(ielr-1)*2+irn
          do 80 ielp=1,nelp
          do 70 ipn=1,3
          ip=(ielp-1)*2+ipn
          iel=(ielr-1)*nelp+ielp

! do loop over grid points
          do 60 ird=0,1
          do 50 ipd=0,1
          id=ird*2+ipd+1
          il=(irn-1)*12+(ipn-1)*4+id

!  boundary  condition at r=rmax
           if (ir.eq.2*nelr+1) go to 40
                     
! boundary condition dpsi/dp = 0 at p=0 and p=pi
           if (ip.eq.1) then
              if (ipd.eq.1) mc(iel,il)=0
              if (ipd.eq.0) mc(iel,il)=(ir-1)*nrow+ird+1

           else if (ip.eq.2*nelp+1) then
              if (ipd.eq.1) mc(iel,il)=0
              if (ipd.eq.0) mc(iel,il)=(ir-1)*nrow+2+(ip-2)*4+ird+1

           else
!  non-boundary grid points
                mc(iel,il)=(ir-1)*nrow+2+(ip-2)*4+id
           endif
           go to 50

40       continue
         
          if (ip.eq.1) then
              if (ipd.eq.1) mc(iel,il)=0
              if (ipd.eq.0) mc(iel,il)=igrid+1

           else if (ip.eq.2*nelp+1) then
              if (ipd.eq.1) mc(iel,il)=0
              if (ipd.eq.0) mc(iel,il)=igrid+1+(ip-2)*2+1
           else
                mc(iel,il)=igrid+1+(ip-2)*2+ipd+1
           endif
           go to 50

50  continue
60  continue
70  continue
80  continue
90  continue
100 continue

           return
           end subroutine conmat
!_______________________________________________________________________________

!  FEM polynomials 5th order 

    subroutine poly(ngaus,root,phi,phid,phidd)

    implicit real*8 (a-h,o-z)
    dimension phi(6,ngaus),phid(6,ngaus),phidd(6,ngaus),root(ngaus)

         do 20 n = 1, ngaus
                x=root(n)
                x0=1.0d0
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

   20    continue
  return
  end subroutine poly

      SUBROUTINE GAUSS (NGAUS,ALPHA,ROOT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALPHA(8),ROOT(8)
! Modfiy for higher quadrature if needed

      if (ngaus.ne.8) stop

      ALPHA(1)=0.362683783378362D0
      ALPHA(2)=0.313706645877887D0
      ALPHA(3)=0.222381034453374D0
      ALPHA(4)=0.101228536290376D0
      ROOT(1)=0.183434642495650D0
      ROOT(2)=0.525532409916329D0
      ROOT(3)=0.796666477413627D0
      ROOT(4)=0.960289856497536D0
 
      IMAX=NGAUS/2
      DO 200 I=1,IMAX
      ALPHA(NGAUS+1-I)=ALPHA(I)
      ROOT(NGAUS+1-I)=-ROOT(I)
200 continue
       RETURN

       end subroutine gauss
!........



 subroutine pfun(i,x,bal)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension fun(6)
      fun(1)=x**2-(1.25)*x**3-(0.5)*x**4+(0.75)*x**5
      fun(2)=(0.5)*x**2-(0.25)*x**3-(0.25)*x**4+(0.25)*x**5
      fun(3)=1.0-2*x**2+x**4
      fun(4)=x-2*x**3+x**5
      fun(5)=x**2+(1.25)*x**3-(0.5)*x**4-(0.75)*x**5
      fun(6)=-(0.25)*x**2-(0.25)*x**3+(0.25)*x**4+(0.25)*x**5

      bal=fun(i)
      
      end subroutine pfun

subroutine psigaia(n,xil,phil,bek,em)
       implicit real*8 (r-z)
       complex*16::em
       complex*16 bek(:)
       em=0
      !n: elementuaren zenbakia
      !ixip: 1-6 polinomioaren indizea
      !iphip: 1-6 polinomioaren indizea
      !bekind: 1-36 koefizienteen indizea
      !xil= xi koordenatu lokala (-1,1)
      !phil= phi koordenatu lokala (-1,1)
      nprim=int(n/int(30))*36
          
      do ixip=1,6

         do iphip=1,6
      ibek=i+j+(5*i-6)+nprim 
      
      call pfun(ixip,xil,s1)
      call pfun(iphip,phil,s2)
      
      em=em+cmplx(bek(ibek))*s1*s2

          end do
     
      end do
      
      end subroutine psigaia


      end program scatter 



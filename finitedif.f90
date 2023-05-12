program diferfin
use tipoak
use mcf_slineales
use mcf_diagonalizacion
complex(kind=dp),parameter::irud=(0.0_dp,1.0_dp)
real(kind=dp),parameter::pi=acos(-1.0_dp)
integer,parameter:: Nxi=30, Nphi=30, N=30
integer,parameter::Nt=8
integer,parameter:: NNN=Nxi*Nphi, NOD=8*Nphi+12, LDIM=3*NOD+1
real(kind=dp),dimension(NNN,NNN)::psi1,psi2
!complex(kind=dp),dimension(Nxi+10,Nphi+10)::psi1prim
!complex(kind=dp),dimension(NNN+10)::psi0primi

integer,dimension(nnn)::indices
real(kind=dp),dimension(nnn)::bprim
real(kind=dp):: det,x,y,em

real(kind=dp),dimension(NNN,1)::psi0,psikalk,ipiv,psiem
real(kind=dp)::hxi,hxi2,hphi,hphi2,ximax,a,batugai,tau,t,u,angelua,xibal
integer::it,i,j,info, nprim, jprim, iprim
real(kind=dp),dimension(NNN)::xi,phi,d,e,xit
ximax=10*pi
hxi=ximax/Nxi
hphi=pi/Nphi
hxi2=hxi*hxi
hphi2=hphi*hphi
psi1=0.0_dp

psi0=0.0_dp
tau=0.5_dp
a=pi
psikalk=0.0_dp

open(unit=16,file="findifgraf.dat",action="write",status="replace")
do i=1,N
 do j=1,N
   xi(j+N*(i-1))=hxi*(j)
   phi(j+N*(i-1))=hphi*i
end do
end do

do i=1,NNN
psikalk(i+N-1,1)=(cos(ximax*cos(pi/2))**2+sin(ximax*cos(pi/2))**2)
end do


do it=1,Nt
   t=tau*it
   nprim=0
   
  !   psikalk(it+N-1,1)=(cos(ximax*cos(pi/2))**2+sin(ximax*cos(pi/2))**2)

  

   do iprim=1,nnn
   
   nprim=nprim+1
    
   jprim=iprim

    
   u=0.0_dp
! if (xi(jprim)<=a) then
           !u=1000000.0_dp
 !  end if
    
   xibal=xi(iprim)
   
   psi0(nprim,1)=psikalk(jprim,1)*(tau*(1-u)-tau/hxi2-tau/hphi2/(xibal**2)+1)
   psi1(iprim,jprim)=1+tau/hxi**2+tau/hphi2/xibal**2
   
 
   
   if (jprim+1>NNN) then
      
  psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim-1,1)*(tau*(1/(2*hxi2)-1/(2*hxi*xibal)))
    psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim-N,1)*(tau/(2*hphi2*xibal**2))


   psi1(iprim,jprim-1)=tau/(2*hxi2)
   psi1(iprim,jprim-N)=tau/(2*hphi2*xibal**2) 
 
     
   else if (jprim+N>NNN) then
        psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim+1,1)*(tau*(1/(2*hxi2)+1/(2*hxi*xibal)))
  psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim-1,1)*(tau*(1/(2*hxi2)-1/(2*hxi*xibal)))
    psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim-N,1)*(tau/(2*hphi2*xibal**2))


   psi1(iprim,jprim+1)=tau/(2*hxi2)
   psi1(iprim,jprim-1)=tau/(2*hxi2)
   psi1(iprim,jprim-N)=tau/(2*hphi2*xibal**2)

  else if (jprim==1) then
           psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim+1,1)*(tau*(1/(2*hxi2)+1/(2*hxi*xibal)))
    psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim+N,1)*(tau/(2*hphi2*xibal**2))


   psi1(iprim,jprim+1)=tau/(2*hxi2)
   psi1(iprim,jprim+N)=tau/(2*hphi2*xibal**2)
  
   else if (jprim-N<1) then
            psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim+1,1)*(tau*(1/(2*hxi2)+1/(2*hxi*xibal)))
  psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim-1,1)*(tau*(1/(2*hxi2)-1/(2*hxi*xibal)))
    psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim+N,1)*(tau/(2*hphi2*xibal**2))


   psi1(iprim,jprim+1)=tau/(2*hxi2)
   psi1(iprim,jprim-1)=tau/(2*hxi2)
   psi1(iprim,jprim+N)=tau/(2*hphi2*xibal**2)
   
  else
   
  psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim+1,1)*(tau*(1/(2*hxi2)+1/(2*hxi*xibal)))
  psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim-1,1)*(tau*(1/(2*hxi2)-1/(2*hxi*xibal)))
    psi0(nprim,1)=psi0(nprim,1)+psikalk(jprim+N,1)*(tau/(2*hphi2*xibal**2))+psikalk(jprim-N,1)*(tau/(2*hphi2*xibal**2))


   psi1(iprim,jprim+1)=tau/(2*hxi2)
   psi1(iprim,jprim-1)=tau/(2*hxi2)
   psi1(iprim,jprim+N)=tau/(2*hphi2*xibal**2)
   psi1(iprim,jprim-N)=tau/(2*hphi2*xibal**2)

   end if
 end do
 
  
   call print_matrix(psi0,15) 
                  call gaussj(psi1,psi0)
                 
   if (it==Nt) then 
           do i=1,NNN
              x=xi(i)*cos(phi(i))
              y=xi(i)*sin(phi(i))
              em=psi0(i,1)**2
              if (xi(i)<a) then
                      write(unit=16, fmt="(3f12.6)"), x, y, 0.0_dp
               else

              write(unit=16,fmt="(3f12.6)"), x, y, em
                end if
           end do
  end if
  
  psikalk=psi0
  psi0=0.0_dp 
 
  end do


  
    
end program diferfin

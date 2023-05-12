module integralak
use tipoak
use polinomioak
use gauss_legendre


real(kind=dp),dimension(6)::pint,pdint,pddint
real(kind=dp),dimension(6,6):: emaitza

pint(1)=integ(p1)
pint(2)=integ(p2)
pint(3)=integ(p3)
pint(4)=integ(p4)
pint(5)=integ(p5)
pint(6)=integ(p6)


pdint(1)=integ(p1d)
pdint(2)=integ(p2d)
pdint(3)=integ(p3d)
pdint(4)=integ(p4d)
pdint(5)=integ(p5d)
pdint(6)=integ(p6d)

pddint(1)=integ(p1dd)
pddint(2)=integ(p2dd)
pddint(3)=integ(p3dd)
pddint(4)=integ(p4dd)
pddint(5)=integ(p5dd)
pddint(6)=integ(p6dd)

bgai(1)      

do i=1,6
      do j=1,6
      emaitza(i,j)=pint(j)*pddint(i)/h+bgai(i)*pint(j)/h
      emaitza(i,j)=emaitza(i,j)+pint(i)*pdint(j)/h+(1-u)*pint(i)*pint(j)


      end do
end do

contains

      !!!!integralak
        
        !bgaiak(x)

        

        function integ(f)
        real(kind=dp)::a,b,integ
        interface 
        use tipoak
        function f(x)
               real(kind=dp),intent(in)::x
               real(kind=dp)::f
        end function f
        end interface

        a=-1.0_dp
        b=1.0_dp
        call qgauss(f,a,b,integ)
        end function pint

       





end integralak

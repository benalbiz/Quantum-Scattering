module polinomioak
use tipoak
use gauss_legendre
public:: p1,p1d,p1dd,p2,p2d,p2dd,p3,p3d,p3dd,p4,p4d,p4dd,p5,p5d,p5dd,p6,p6d,p6dd
public:: bgaia1,bgaia2,bgaia3,bgaia4,bgaia5,bgaia6
public:: derb1,derb2,derb3,derb4,derb5,derb6
contains

      !!!!!integratu funtzio guztiak. ez ahaztu bgaia be definitzie. ta gero if batzukin juan batzen emaitzak kasu bakoitzerako.

      function p1(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p1
      p1=x**2-(1.25_dp)*x**3-(0.5_dp)*x**4+(0.75_dp)*x**5
      end function p1
      function p2(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p2
      p2=(0.5_dp)*x**2-(0.25_dp)*x**3-(0.25_dp)*x**4+(0.25_dp)*x**5
      end function p2

      function p3(x)
       real(kind=dp),intent(in)::x
      real(kind=dp)::p3
      p3=1.0_dp-2*x**2+x**4
      end function p3

      function p4(x)
       real(kind=dp),intent(in)::x
       real(kind=dp)::p4
      p4=x-2*x**3+x**5
      end function p4

      function p5(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p5
      p5=x**2+(1.25_dp)*x**3-(0.5_dp)*x**4-(0.75_dp)*x**5
      end function p5

      function p6(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p6
      p6=-(0.25_dp)*x**2-(0.25_dp)*x**3+(0.25_dp)*x**4+(0.25_dp)*x**5
      end function p6

      function p1d(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p1d
      p1d=2*x-(3.75_dp)*x**2-2*x**3+(3.75_dp)*x**4
      end function p1d

      function bgaia1(xi0,hxi,x)
      real(kind=dp),intent(in)::x,xi0,hxi
      real(kind=dp)::bgaia1
      bgaia1=1/(xi0+hxi*x)*(2*x-(3.75_dp)*x**2-2*x**3+(3.75_dp)*x**4)/hxi
      end function bgaia1


      function p2d(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p2d
      p2d=(0.5_dp)*x-(0.75_dp)*x**2-x**3+(1.25_dp)*x**4
      end function p2d

      function bgaia2(xi0,hxi,x)
      real(kind=dp),intent(in)::x,xi0,hxi
      real(kind=dp)::bgaia2
      bgaia2=1/(xi0+hxi*x)*((0.5_dp)*x-(0.75_dp)*x**2-x**3+(1.25_dp)*x**4)/hxi
      end function bgaia2

      function p3d(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p3d
      p3d=-4*x+4*x**3 
      end function p3d

      function bgaia3(xi0,hxi,x)
      real(kind=dp),intent(in)::x,xi0,hxi
      real(kind=dp)::bgaia3
      bgaia3=1/(xi0+hxi*x)*(-4*x+4*x**3)/hxi
      end function bgaia3
 

      function p4d(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p4d
      p4d=1.0_dp-6*x**2+5*x**4
      end function p4d

      function bgaia4(xi0,hxi,x)
      real(kind=dp),intent(in)::x,xi0,hxi
      real(kind=dp)::bgaia4
      bgaia4=1/(xi0+hxi*x)*(1.0_dp-6*x**2+5*x**4)/hxi
      end function bgaia4

      function p5d(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p5d
      p5d=2*x+(3.75_dp)*x**2-2*x**3-(3.75_dp)*x**4
      end function p5d

      function bgaia5(xi0,hxi,x)
      real(kind=dp),intent(in)::x,xi0,hxi
      real(kind=dp)::bgaia5
      bgaia5=1/(xi0+hxi*x)*(2*x+(3.75_dp)*x**2-2*x**3-(3.75_dp)*x**4)/hxi
      end function bgaia5

      function p6d(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::p6d
      p6d=-(0.5_dp)*x-(0.75_dp)*x**2+x**3+(1.25_dp)*x**4
      end function p6d
      
      function bgaia6(xi0,hxi,x)
      real(kind=dp),intent(in)::x,xi0,hxi
      real(kind=dp)::bgaia6
      bgaia6=1/(xi0+hxi*x)*(-(0.5_dp)*x-(0.75_dp)*x**2+x**3+(1.25_dp)*x**4)/hxi
      end function bgaia6

      function p1dd(xi0,hxi,hphi,x)
      real(kind=dp),intent(in)::x,xi0,hxi,hphi
      real(kind=dp)::p1dd
      p1dd=(2.0_dp-(7.5_dp)*x-6*x**2+15*x**3)/(hphi**2*(xi0+hxi*x)**2)
      end function p1dd

      function derb1(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::derb1
      derb1=(2.0_dp-(7.5_dp)*x-6*x**2+15*x**3)
      end function derb1


      function p2dd(xi0,hxi,hphi,x)
      real(kind=dp),intent(in)::x,xi0,hxi,hphi
      real(kind=dp)::p2dd
      p2dd=((1.0_dp/2.0_dp)-(1.5_dp)*x-3*x**2+5*x**3)/(hphi**2*(xi0+hxi*x)**2)

      end function p2dd

      function derb2(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::derb2
      derb2=((1.0_dp/2.0_dp)-(1.5_dp)*x-3*x**2+5*x**3)

      end function derb2


      function p3dd(xi0,hxi,hphi,x)
      real(kind=dp),intent(in)::x,xi0,hxi,hphi
      real(kind=dp)::p3dd
      p3dd= (-4.0_dp+12*x**2)/(hphi**2*(xi0+hxi*x)**2)
      end function p3dd

      function derb3(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::derb3
      derb3= (-4.0_dp+12*x**2)
      end function derb3


      function p4dd(xi0,hxi,hphi,x)
      real(kind=dp),intent(in)::x,xi0,hxi,hphi
      real(kind=dp)::p4dd
      p4dd=(-12*x+20*x**3)/(hphi**2*(xi0+hxi*x)**2)
      end function p4dd

      function derb4(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::derb4
      derb4=(-12*x+20*x**3)
      end function derb4


      function p5dd(xi0,hxi,hphi,x)
      real(kind=dp),intent(in)::x,xi0,hxi,hphi
      real(kind=dp)::p5dd
      p5dd=(2.0_dp+(7.5_dp)*x-6*x**2-15*x**3)/(hphi**2*(xi0+hxi*x)**2)
      end function p5dd

      function derb5(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::derb5
      derb5=(2.0_dp+(7.5_dp)*x-6*x**2-15*x**3)
      end function derb5



      function p6dd(xi0,hxi,hphi,x)
      real(kind=dp),intent(in)::x,xi0,hxi,hphi
      real(kind=dp)::p6dd
      p6dd=(-(0.5_dp)-(1.5_dp)*x+3*x**2+5*x**3)/(hphi**2*(xi0+hxi*x)**2)

      end function p6dd

      function derb6(x)
      real(kind=dp),intent(in)::x
      real(kind=dp)::derb6
      derb6=(-(0.5_dp)-(1.5_dp)*x+3*x**2+5*x**3)

      end function derb6

end module polinomioak 

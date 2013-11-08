subroutine poisson_solve(rho,phi)
  implicit none
  include 'param.h'
  double precision :: rho(ns), phi(ns)             !!!Change
  double precision :: a, G, K, gamma, deltac,xl,xr, C
  double precision, allocatable :: x(:),dphi(:)
  integer :: n,i,test, count
  double precision :: h  


  
  n=ns       !!!Change

  xl=0d0
  xr=1
  
!!!!Stepsize    
  h=(xr-xl)/(n-1)!!!!!!Include boundery points in array


!!!!Allocate arrays  

  allocate (x(n))
  allocate (dphi(n))

 
  do i=1,n
    dphi(i)=0
  end do
  
  
  do i=1,n
    x(i)=(i-1)*h
  end do


  !!!!First Poisson solved, trial phi
  phi(1)=0
  phi(2)=0
  phi(n)=0


  do i=1,n-1
    phi(n)=phi(n)+h*rho(i)/((n-i)*h)
  end do
  
!  print*,phi(n)
 
  do i=3,n-1
!    phi(i+1)=(4*3.14*h**2*x(i)**2*rho(i)+(x(i)-h/2)**2*(phi(i)-phi(i-1))+(x(i)+h/2)**2*phi(i))/(x(i)+h/2)**2
     phi(i)= phi(n)/(1-h)*x(i)
  end do
  

  
  do i=1,n
    dphi(i)=phi(i)
  end do
  
  !!!!Poisson Relaxation
  test=1
  count=0
  do while (test==1)
    
    do i=2,n-1
      phi(i)=(phi(i+1)*(x(i)+h/2)**2+phi(i-1)*(x(i)-h/2)**2-4*3.14*h**2*x(i)**2*rho(i))/((x(i)+h/2)**2+(x(i)-h/2)**2)
    end do
    
    phi(1)=phi(2)
    
    test=0
    do i=1,n

      dphi(i)=abs(dphi(i)-phi(i))
      if (dphi(i).gt.1d-10) then
        test=1
      end if 
    end do  

    do i=1,n
      dphi(i)=phi(i)
    end do  
  

  count =count +1
  end do
!  print*, "Internal Iterations",count
  
!  print*,"phi"
!  do i=1,n  
!    print*, phi(i)
!  end do 
!  print*,count

end subroutine poisson_solve

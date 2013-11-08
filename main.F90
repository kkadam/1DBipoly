program main
  implicit none
  include 'param.h'
  double precision :: a, G, Rho_0, gamma1, gamma2, rc, dR, K1, K2, C, Phi_s, Phi_c , C1
  double precision :: Rho1_c, Rho2_c, deltaC, prevC, deltaK1,k1old, deltaK2, k2old
  double precision, allocatable :: R1(:), R2(:), R(:), Rho1(:), Rho2(:), Rho(:)
  double precision, allocatable :: Phi(:), x(:), H1(:), H2(:), H(:)
  integer :: i , conv,counter,n1,n2
  double precision :: cpu1, cpu2, cput
  
  print*,"===========================================================================" 
  print*, "SCF started!"

  call cpu_time(cpu1)  
  
  a=1
  G=1
  Rho_0=1

  gamma1=1+1.0/np1
  gamma2=1+1.0/np2
  
  n1=nc
  n2=ns-nc+1
  
  rc=(Nc-1)/(Ns-1)*a
  dR=a/(Ns-1)
  
  
  allocate (R1(N1))
  allocate (R2(N2))
  allocate (R(Ns))
  allocate (Rho1(N1))
  allocate (Rho2(N2))
  allocate (Rho(Ns))
  allocate (x(Ns))
  allocate (Phi(Ns))
  allocate (H(Ns))
  
!!!!Get Position Array  
  do i=1, Ns
    x(i)=(i-1)*dR
  end do


	
!!!!Assume rho  
  do i=1,Ns
    Rho(i)= (1-x(i))
  end do  
  
  do i=1, N1
    Rho1(i)= Rho(i)
  enddo
  
  do i=1, N2
    Rho2(i)=Rho(N1+i-1)
  enddo
     
!  call printdefault(rho)
!  call findpar(rho,2.27, 0.30)
  counter=1
  deltaK1=1.0
  deltaK2=1.0 
  K1old=0.0
  K2old=0.0
!!!!Iterate Till Convergence!!!!!!!!!  
  do while ((deltaK1.gt.1d-4).and.(deltaK2.gt.1d-4))
!      call mass_solve(rho,phi)
      call poisson_solve(rho,phi)
      
      phi_c=phi(N1)
      phi_s=phi(Ns)
      
      rho1_c=rho1(N1)
      rho2_c=rho2(1)
      
     rho2_c= rho1_c*mu2/mu1

       K2=(phi_s-phi_c)/(np2+1)/rho2_c**(1.0/np2)
       K1=K2*rho2_c**(gamma2)/rho1_c**(gamma1)
      
       C1 = phi_c+(np1+1)*K2*rho2_c**gamma2/rho1_c
      
!       if (counter==1) then
!          deltaK1=1
!          deltaK2=1
!       else
      	  deltaK1=Abs(K1-K1old)
	  deltaK2=Abs(K2-K2old)
!       endif  
      
       K1old=K1
       K2old=K2
      
      do i=1,N2
        Rho2(i)=(Phi_s/K2/(np2+1))**np2*(1-Phi(N1+i-1)/Phi_s)**np2
      enddo
     
!     print*, "const2",(Phi_s/K2/(np2+1))**np2
     
      do i=1,N1
        Rho1(i)=(C1/K1/(np1+1))**np1*(1-Phi(i)/C1)**np1
      enddo
     
!     print*, "const1",(C1/K1/(np1+1))**np1

     
      Rho_0=Rho1(1)
     
      do i=1,N1
        Rho(i)=Rho1(i)/Rho_0
      enddo
  
      do i=1,N2-1
        Rho(N1+i)=Rho2(i+1)/Rho_0
      enddo
    
      Print*, counter, "th iteration"
      print*, "K1 = ",K1, "K2 = ", K2
      print*, "deltaK1 = ", deltaK1, "deltaK2 = ", deltaK2
      counter=counter+1
  enddo
     
  call cpu_time(cpu2)
  cput=(cpu2-cpu1)/60.0
  

  call findpar(rho,phi,K1,K2,counter,cput)
 
end program main  
  


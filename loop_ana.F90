program ana
!Analytical solution for n1=5 and n2=1 bipolytrope
  implicit none
  double precision :: a, G, Rho_0, gamma1, gamma2, rc, dR, K1, K2, C, Phi_s , C1, mu1, mu2, xi_c
  double precision :: Rho1_c, Rho2_c, deltaC, prevC, np1, np2, dXi
  double precision :: theta_c,gradTheta_c,phi_c,eta_c,gradPhi_c,Lambda_c,B,eta_s,Pi, dEta
  double precision, allocatable :: R1(:), R2(:), R(:), Rho1(:), Rho2(:), Rho(:),eta(:)
  double precision, allocatable :: Phi2(:), x(:), xi(:)
  integer :: Ns, Nc,Ne, i ,n1,n2, conv,counter,j
  double precision :: Mcore, Mt, Mc, gradPhi_s
  

  
!Specify the two atomic weights and core radius 
  mu2=0.75
  mu1=1
  

  
!Specify other constants and gridsize for the core  
  Pi=3.14159265359
  G=1
  Rho_0=1
  K1=1.0
  Nc=101     
  
  
!Allocate arrays for core  
  allocate (Rho1(Nc))
  allocate (Xi(Nc))
  allocate (R1(Nc))

!Specify envelope resolution  
  Ne=101


!Allocate envelope arrays  
  allocate (R2(Ne))
  allocate (Rho2(Ne))
  allocate (Phi2(Ne))
  allocate(eta(Ne))
  
do j=1,100
  Xi_c=0.1*j  
  
!Get Xi array  
  dXi=Xi_c/(Nc-1)
 
  do i=1, Nc
    Xi(i)=(i-1)*dXi
  enddo

!Get core radius array  
  do i=1,Nc
    R1(i)= (K1/(G*rho_0**(4.0/5)))**0.5*(3/(2*Pi))**0.5*Xi(i)
  enddo

  
!Get core density profile  
  do i=1, Nc
    Rho1(i)= Rho_0*(1+Xi(i)**2/3.0)**(-2.5)
  enddo
  
!Constants from interface conditions 
  theta_c=(1+Xi_c**2/3)**(-0.5)
  gradTheta_c=-Xi_c/3*(1+Xi_c**2/3)**(-1.5)  
  
  phi_c=1.0
  eta_c=3**(0.5)*xi_c*(mu2/mu1)*theta_c**2  
  gradPhi_c=3**(0.5)*theta_c**(-3.0)*gradTheta_c   
  Lambda_c=1/eta_c+1/phi_c*gradPhi_c
  
  
  B= eta_c-Pi/2+atan(Lambda_c)   
  eta_s=Pi+B
  A=phi_c*eta_c*(1+Lambda_c**2.0)**(0.5)
  
  Rho1_c=Rho1(Nc)
  Rho2_c=(mu2/mu1)*Rho1_c
  
  dEta=(eta_s-eta_c)/(Ne)

!Get eta array  
  eta(1)=eta_c
  do i=2, Ne
    eta(i)=eta_c+dEta*i 
  enddo

!Get envelope density array		
  do i=1, Ne
    Phi2(i)=A*sin(eta(i)-B)/eta(i)
  enddo  
  
  Rho2(1)=Rho2_c
  do i=2, Ne
    Rho2(i)=Rho_0*(mu2/mu1)*Theta_c**5*Phi2(i)
  enddo
  
  
  do i=1, Ne
    R2(i)=(1.0)*(mu1/mu2)/(Theta_c**2)/(2*Pi)**(0.5)*eta(i)
  enddo
  
!   print*, "theta_c", theta_c
!  print*, "gradTheta_c", gradTheta_c
!   print*, "phi_c",phi_c

! print*,"Lambda_c",Lambda_c
! print*,"A",A
! print*,"B",B
!  print*, "eta_c",eta_c
!  print*, "eta_s",eta_s
!    print*, "gradphi_c",gradphi_c

  
  
!Find core mas ratio
  Mcore= (K1**3/G**3/Rho_0**(2.0/5)*6/Pi)**(0.5)*(Xi_c**3*(1+Xi_c**2/3.0)**(-1.5))
  
  gradPhi_s=A/eta_s**2*(eta_s*cos(eta_s-B)-sin(eta_s-B))
  Mt= (K1**3/G**3/Rho_0**(2.0/5)*2/Pi)**(0.5)*(mu1/mu2)**2.0/theta_c*(-eta_s**2*gradPhi_s)
  
  Mc=Mcore/Mt
  Rc=R1(Nc)/R2(Ne)
  
!  print*, "gradPhi_s",gradPhi_s 
!  print*, "Mcore", Mcore
!  print*, "Mt", Mt
  print*, "Mc=", Mc,"Rc=",Rc
  
  
  
!  open(unit=14,file='test.dat')	 !!!Change
!    do i=0,Nc
!      write(14,*) R1(i), Rho1(i)
!    enddo
!    do i=2,Ne
!      write(14,*) R2(i), Rho2(i)
!    enddo
!   close(14)
!  print*, "File test printed."      !!!Change

    open(unit=12,file='SC.dat',access='APPEND')  
    write(12,*) Rc,Mc
    close(12)
	
enddo
	
	
end program ana  

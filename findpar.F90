subroutine findpar(rho,pot, K1, K2, count, cput, ns, nc)
  implicit none
  include 'param.h'
 

   double precision, dimension(ns) :: rho, press, pot
  
   double precision :: m, dr,pi,r,rho_2i,m_core,q,vol,presse, K1, K2, W, rcore,&
   pmax
   integer :: i,j, count, ns ,nc
  
  character(len=100) :: img_file, info_file
  character*50 char_np1, char_nc, char_ns, char_np2, char_mu1, char_mu2, &
  char_rcore, char_mcore, char_m, char_q, char_vol, char_W, char_3P,&
  char_pmax, char_count, char_cput, char_K1,char_K2

  double precision :: cput
   
   Pi=3.14159265359 
   dr=1.0/(ns-1)
   
   rcore=(nc-0.5)/(ns-0.5)
   
   rho_2i=rho(nc+1)

!Find mass  
   m=0.0
   m_core=0.0
   
   do i=1,ns   
      r=(i-0.5)*dr        
      m=m+4*pi*r**2*rho(i)*dr
      if (rho(i).gt.rho_2i) then         
        m_core=m_core+4*pi*r**2*rho(i)*dr
      endif
   enddo

   q=m_core/m
   
   
!Find volume  
   vol=0.0
   do i=1,ns 
     r=(i-1.5)*dr
     vol=vol+4*pi*r**2*dr
   enddo
!Find potential energy W
	
  W=0.0
  do i=1,ns
    r=(i-0.5)*dr
    W=W-0.5*pot(i)*4*pi*r**2*rho(i)*dr
  enddo   	     

  
!Find pressure energy P
  do i=1,nc
    press(i)=K1*rho(i)**(1+1.0/np1)
  enddo

  do i=nc+1,ns
    press(i)=K2*rho(i)**(1+1.0/np2)
  enddo
	
  pressE=0.0
  do i=1,ns
    r=(i-0.5)*dr
    pressE=pressE+press(i)*4*pi*r**2*dr
  enddo    

  pmax=maxval(press)
  
!Write files
!Convert to strings
  if (nc.lt.100) then
    write (char_nc, "(I2)") nc
  else
    write (char_nc, "(I3)") nc
  endif

  if ((ns.lt.1000).and.(ns.gt.100)) then
    write (char_ns, "(I3)") ns
  else
    write (char_ns, "(I4)") ns
  endif 

  write (char_np1, "(F3.1)") np1
  write (char_np2, "(F3.1)") np2
  write (char_mu1,"(F3.1)") mu1
  write (char_mu2,"(F3.1)") mu2

  write (char_rcore,"(F6.4)") rcore
  write (char_mcore,"(F6.4)") m_core
  write (char_m,"(F6.4)") m
  write (char_q,"(F6.4)") q
  write (char_vol,"(F6.4)") vol
  write (char_W,"(F6.4)") W
  write (char_3P,"(F6.4)") presse*3.0
  write (char_pmax,"(F6.4)") pmax 
  write (char_count,"(I3)") count
  write (char_cput,"(F6.4)") cput
  write (char_rcore,"(F6.4)") rcore 
  write (char_K1,"(F6.4)") K2 
  write (char_K2,"(F6.4)") K2 

!Write rho file  
  img_file='Bi_'//trim(char_nc)//"_"//trim(char_ns)//"_"//trim(char_np1)//'w'&
  //trim(char_mu1)//'_'//trim(char_np2)//'w'//trim(char_mu2)   

  open(unit=10,file=img_file)
    do i=1,ns  
      write(10,*) rho(i) 
    enddo
  close(10)   
  
!Write info files  
  info_file='Bi_'//trim(char_nc)//"_"//trim(char_ns)//"_"//trim(char_np1)//'w'&
  //trim(char_mu1)//'_'//trim(char_np2)//'w'//trim(char_mu2)//".info"

  
  open(unit=10,file=info_file)
  write(10,*) trim(char_np1)," ",trim(char_np2)," ",trim(char_mu1)," ",trim(char_mu2),&
  " ",trim(char_nc)," ",trim(char_ns)," ",trim(char_rcore)," ",trim(char_mcore),&
  " ",trim(char_m)," ",trim(char_q)," ",trim(char_vol)," ",trim(char_W)," ",&
  trim(char_3P)," ",trim(char_pmax)," ",trim(char_count)," ",trim(char_cput)
  close(10)  
  

  open(unit=12,file='SC_HR.dat',access='APPEND')  
    write(12,*) trim(char_rcore)," ", trim(char_q)
  close(12)
  
  
  print*,"================================SUMMARY===================================="
  print*,"n_core  = ", trim(char_np1), "  n_env = ", trim(char_np2)
  print*,"mu_core = ", trim(char_mu1), "  mu_env = ", trim(char_mu2)
  print*, "Core mass = ", trim(char_mcore)
  print*, "Core mass fraction = ", trim(char_q)
  print*,"Resolution = ", trim(char_ns)
  print*,"Core Resolution = ", trim(char_nc)  
  print*,"  rcore  ","  M     ", "  V   ","    -W   ","  3PI  ","   P_max  ","   K1   "&
  ,"   K2   "
  print*,trim(char_rcore),"   ",trim(char_m),"  ", trim(char_vol),"  ",trim(char_W)&
  ,"  ",trim(char_3P),"  ",trim(char_pmax),"  ",trim(char_K1),"  ",trim(char_K2)
  
  print*,"cpu time =", trim(char_cput) , " min"
  
  print*,"==============================OUTPUT FILES================================="    
  
  print*,trim(info_file)
  print*,trim(img_file)   
  print*,"===========================================================================" 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
!   print*,"findmass mass=", m, "masscount", count 
end subroutine findpar
   

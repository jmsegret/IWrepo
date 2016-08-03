!-builds subd.dat files 

      Program buildsubd
!********************************


      IMPLICIT NONE

      
      INTEGER :: nsubd      !number of subdomains
      REAL :: Htotal     !height of domain
      REAL :: k,i
      INTEGER :: j

      write(*,*) 'Enter the number of subdomains'
      read(*,*) nsubd

      open(66,file='subd.dat')
      write(66,"(I2)") nsubd

      k=0.0
      Do j =1,nsubd
      write(*,*) 'Enter height of the next subdomain (bottom to top)'
      read(*,*) i
      write(66,"(F5.3,' ',F5.3)") 0.0+k, i
      k=k+i
      write(*,*) 'The current height of your subdomain is' , k
      write(*,'(A,I2,A)')'You have ',nsubd-j, ' subdomains remaining'
      end do


      End Program buildsubd
 
      

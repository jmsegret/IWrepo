!----------------------------------------------------------
!> A module to generate pseudo-random numbers.
!>
!> \author Kris R.
!> \date June 23, 2016
!>
!----------------------------------------------------------
      module rng

         use iso_fortran_env
         implicit none
         save

         !> The initial seed
         integer(int64) :: seed = 0

         !> The current integer value of the rng
         integer(int64) :: current = 0
         
      contains

         !----------------------------------------
         !> Seed the rng by xoring the time and procid.
         !> This must be called before using the rng
         !> or else an error will be thrown.
         !---------------------------------------
         subroutine seed_generator()
            save         

            integer(int64) :: mytime, myprocid

            call system_clock(mytime)
            myprocid = getpid()
            seed = ieor(mytime,myprocid)

            current = seed

            return

         end subroutine seed_generator
         
         !----------------------------------------
         !> This implements the xorshift rng.
         !> Better generators are available, however
         !> this will be sufficient for our work and
         !> is better than the lcg generator used
         !> previously. Additionally, this means we
         !> don't rely on the native FORTRAN rand()
         !> function which has different implementations
         !> on different compilers.
         !>
         !> The parameters in this function are chosen
         !> careful so do not change them unless you
         !> really know what you're doing.
         !>
         !> \return A uniformly distributed (pseudo) random
         !> real number in the interval (0,1)
         !>
         !> \throws seed_generator() must be called before using
         !> this function or else an error will be thrown.
         !---------------------------------------
         function random()
            save

            real(real64) :: random
            integer(int64), parameter :: a1=21, a2=-35, a3=4
            real(real64), parameter :: q=5.42101086242752217e-20

            if(seed.eq.0) then
               print *, "RNG has not been seeded!"
               stop
            endif
  
            !perform xorshift
            current = ieor(current, ishft(current, a1))
            current = ieor(current, ishft(current, a2))
            current = ieor(current, ishft(current, a3))
            random = current*q + 0.5

         end function random
         
         
         !----------------------------------------
         !> Return the value used to seed the rng.
         !---------------------------------------
         function get_seed()
            
            integer(int64) :: get_seed

            get_seed = seed

         end function get_seed

         

      end module rng
!----------------------------------------------------------

!----------------------------------------------------------
!- Test harness for this module.
#if 0
      program test_random

         use rng
         implicit none
         
         integer, parameter :: N = 50000
         integer :: i
         integer(int64) :: myseed
         real(real64) :: myrandom(3) = 0.0

         myseed = get_seed()
         print *, "myseed: ", myseed

         call seed_generator()
         myseed = get_seed()
         print *, "myseed: ", myseed

         print*, "myrandom: ", myrandom
         open(77,file='test_random.dat')
         do i=1,N
            
            myrandom(1) = random()
            myrandom(2) = random()
            myrandom(3) = random()
            write(77,*) myrandom

         enddo !i=1,N
         close(77)

      end program test_random

#endif
!----------------------------------------------------------

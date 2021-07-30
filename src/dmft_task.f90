!!========================================================================
!!>>> driver subroutines: layer 1                                      <<<
!!========================================================================

!!
!! @sub dmft_driver
!!
!! core subroutine, dispatch all the computational tasks
!!
  subroutine dmft_driver()
     use control, only : task

     implicit none

! we have to preprocess the self-energy functions at first
! calculate: sigma -> sigoo -> sigoo - sigdc
     call cal_sigoo()

! calculate: sigma -> sigma - sigdc (sigma is updated)
     call cal_sigma()

! main scheduler
! task = 1, calculate hybridization function, for one-shot calculation
! task = 2, calculate density correction, for self-consistent calculation
! task = 3, calculate fermi level
! task = 4, calculate impurity levels
! task = 5, calculate complex dft + dmft eigenvalues
! task = 6, calculate lattice green's functions
! task = 999, only for test
     DISPATCHER: select case ( task )
         !
         case (1)
             call dmft_try1()
         !
         case (2)
             call dmft_try2()
         !
         case (3)
             call dmft_try3()
         !
         case (4)
             call dmft_try4()
         !
         case (5)
             call dmft_try5()
         !
         case (6)
             call dmft_try6()
         !
         case (999)
             call dmft_try999()
         !
         case default
             call s_print_error('dmft_driver','this feature is not supported')
         !
     end select DISPATCHER

     return
  end subroutine dmft_driver

!!========================================================================
!!>>> driver subroutines: layer 2                                      <<<
!!========================================================================

!!
!! @sub dmft_try1
!!
!! to calculate the local green's function, generate key inputs for the
!! quantum impurity solvers. the fermi level may be updated, depending
!! on the configuration parameter. this subroutine is suitable for the
!! one-shot dft + dmft calculations.
!!
  subroutine dmft_try1()
     use constants, only : dp, zero
     use constants, only : mystd

     use control, only : cname
     use control, only : lfermi
     use control, only : fermi
     use control, only : myid, master

     use context, only : eimps, eimpx
     use context, only : green
     use context, only : weiss, delta

     implicit none

! local variables
! lattice occupancy
     real(dp) :: occup

! try to search the fermi level
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Fermi'
     endif ! back if ( myid == master ) block
     !
     occup = zero
     !
     if ( lfermi .eqv. .true. ) then
         call cal_fermi(occup)
     else
         if ( myid == master ) then
             write(mystd,'(4X,a)') 'SKIP'
         endif ! back if ( myid == master ) block
     endif ! back if ( lfermi .eqv. .true. ) block
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the local impurity levels
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Level'
     endif ! back if ( myid == master ) block
     !
     call cal_eimps()
     !
     call cal_eimpx()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the local green's function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Green'
     endif ! back if ( myid == master ) block
     !
     call cal_green()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the local weiss's function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Weiss'
     endif ! back if ( myid == master ) block
     !
     call cal_weiss()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the hybridization function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Hybri'
     endif ! back if ( myid == master ) block
     !
     call cal_delta()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save fermi...'
         call dmft_dump_fermi(fermi, occup, zero)
         !
         write(mystd,'(4X,a)') 'save eimps...'
         call dmft_dump_eimps(eimps)
         !
         write(mystd,'(4X,a)') 'save eimpx...'
         call dmft_dump_eimpx(eimpx)
         !
         write(mystd,'(4X,a)') 'save green...'
         call dmft_dump_green(green)
         !
         write(mystd,'(4X,a)') 'save weiss...'
         call dmft_dump_weiss(weiss)
         !
         write(mystd,'(4X,a)') 'save delta...'
         call dmft_dump_delta(delta)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try1

!!
!! @sub dmft_try2
!!
!! to calculate the density correction, generate key inputs for the dft
!! engine. the fermi level may be updated, depending on the configuration
!! parameter. this subroutine is suitable for the fully self-consistent
!! dft + dmft calculations.
!!
  subroutine dmft_try2()
     use constants, only : dp, zero
     use constants, only : mystd

     use control, only : cname
     use control, only : lfermi
     use control, only : fermi
     use control, only : myid, master

     use context, only : gamma

     implicit none

! local variables
! lattice occupancy
     real(dp) :: occup

! correction to band energy
     real(dp) :: ecorr

! try to search the fermi level
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Fermi'
     endif ! back if ( myid == master ) block
     !
     occup = zero
     !
     if ( lfermi .eqv. .true. ) then
         call cal_fermi(occup)
     else
         if ( myid == master ) then
             write(mystd,'(4X,a)') 'SKIP'
         endif ! back if ( myid == master ) block
     endif ! back if ( lfermi .eqv. .true. ) block
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the density correction
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Gamma'
     endif ! back if ( myid == master ) block
     !
     ecorr = zero
     !
     call cal_gamma(ecorr)
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save fermi...'
         call dmft_dump_fermi(fermi, occup, ecorr)
         !
         write(mystd,'(4X,a)') 'save gamma...'
         call dmft_dump_gamma(gamma)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try2

!!
!! @sub dmft_try3
!!
!! try to search the fermi level, the global variable `fermi` may be
!! updated in this subroutine. it is just for testing purpose.
!!
  subroutine dmft_try3()
     use constants, only : dp, zero
     use constants, only : mystd

     use control, only : cname
     use control, only : lfermi
     use control, only : fermi
     use control, only : myid, master

     implicit none

! local variables
! lattice occupancy
     real(dp) :: occup

! check lfermi at first
     call s_assert2(lfermi .eqv. .true., 'lfermi must be true')

! try to search the fermi level
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Fermi'
     endif ! back if ( myid == master ) block
     !
     occup = zero
     !
     call cal_fermi(occup)
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save fermi...'
         call dmft_dump_fermi(fermi, occup, zero)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try3

!!
!! @sub dmft_try4
!!
!! try to determine the local impurity levels. this subroutine is just
!! for testing purpose.
!!
  subroutine dmft_try4()
     use constants, only : mystd

     use control, only : cname
     use control, only : myid, master

     use context, only : eimps, eimpx

     implicit none

! try to calculate the local impurity levels
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Level'
     endif ! back if ( myid == master ) block
     !
     call cal_eimps()
     !
     call cal_eimpx()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save eimps...'
         call dmft_dump_eimps(eimps)
         !
         write(mystd,'(4X,a)') 'save eimpx...'
         call dmft_dump_eimpx(eimpx)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try4

!!
!! @sub dmft_try5
!!
!! try to calculate all the complex dft + dmft eigenvalues. the subroutine
!! can be used in the postprocessing procedure.
!!
  subroutine dmft_try5()
     use constants, only : dp, mystd

     use control, only : cname
     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : myid, master

     use context, only : qbnd

     implicit none

! local variables
! status flag
     integer  :: istat

! dummy array, used to save the eigenvalues of H + \Sigma(i\omega_n)
     complex(dp), allocatable :: eigs(:,:,:,:)

! dummy array, used to save the eigenvalues of H + \Sigma(ioo)
     complex(dp), allocatable :: einf(:,:,:)

! allocate memory
     allocate(eigs(qbnd,nmesh,nkpt,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('dmft_try5','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(einf(qbnd,nkpt,nspin),       stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('dmft_try5','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! try to diagonalize the effective hamiltonian
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Eigen'
     endif ! back if ( myid == master ) block
     !
     call cal_eigsys(eigs, einf)
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save eigen...'
         call dmft_dump_eigen(eigs)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     if ( allocated(eigs) ) deallocate(eigs)
     if ( allocated(einf) ) deallocate(einf)

     return
  end subroutine dmft_try5

!!
!! @sub dmft_try6
!!
!! try to calculate all the physical quantities related to momentum, such
!! as lattice green's functions, and momentum-resolved spectral functions.
!! the subroutine can be used in the postprocessing procedure.
!!
  subroutine dmft_try6()
     implicit none

     return
  end subroutine dmft_try6

!!
!! @sub dmft_try999
!!
  subroutine dmft_try999()
     implicit none

     return
  end subroutine dmft_try999

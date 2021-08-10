!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : cal_sigoo
!!!           cal_sigma
!!!           cal_fermi
!!!           cal_eimps
!!!           cal_eimpx
!!!           cal_green
!!!           cal_green_tetra
!!!           cal_weiss
!!!           cal_delta
!!!           cal_gamma
!!! source  : dmft_flow.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/29/2021 by li huang (created)
!!!           08/09/2021 by li huang (last modified)
!!! purpose : implement the main work flow of dft + dmft calculation.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub cal_sigoo
!!
!! try to calculate the asymptotic values for self-energy functions. then
!! the double-counting terms will be removed as well. this function works
!! for Matsubara self-energy functions (bare) only.
!!
  subroutine cal_sigoo()
     use constants, only : dp
     use constants, only : czero

     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh

     use context, only : qdim
     use context, only : sigdc, sigoo, sigma

     implicit none

!! local parameters
     ! how many frequency points are included to calculate the asymptotic
     ! values of Matsubara self-energy function.
     integer, parameter :: mcut = 16

!! local variables
     ! loop index for impurity sites
     integer :: t

     ! loop index for spins
     integer :: s

     ! loop index for frequency mesh
     integer :: m

     ! status flag
     integer :: istat

     ! dummy array for the indices of frequency points
     integer, allocatable :: ip(:)

     ! dummy array for the Matsubara self-energy functions
     complex(dp), allocatable :: Sm(:,:)

!! [body

     ! allocate memory
     allocate(ip(mcut), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sigoo','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Sm(qdim,qdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sigoo','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! check working axis.
     ! now only the Matsubara frequency axis is supported.
     call s_assert2(axis == 1, 'axis is wrong')

     ! reset sigoo
     sigoo = czero

     ! build integer array for indices of frequency points
     call s_linspace_i(nmesh + 1 - mcut, nmesh, mcut, ip)

     ! loop over quantum impurities, spins, and frequency points.
     !
     ! we count the last `mcut` frequency points, then we try to
     ! calculate the averaged values. up to now, the double counting
     ! terms have not been substracted from sigma. in other words,
     ! sigma is still bare.
     do t=1,nsite
         do s=1,nspin
             Sm = czero
             !
             do m=1,mcut
                 Sm = Sm + sigma(:,:,ip(m),s,t)
             enddo ! over m={1,mcut} loop
             !
             sigoo(:,:,s,t) = Sm / float(mcut)
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

     ! we substract the double counting terms from sigoo
     sigoo = sigoo - sigdc

     ! deallocate memory
     if ( allocated(ip) ) deallocate(ip)
     if ( allocated(Sm) ) deallocate(Sm)

!! body]

     return
  end subroutine cal_sigoo

!!
!! @sub cal_sigma
!!
!! try to substract the double counting terms from the bare Matsubara
!! self-energy functions. this function works for Matsubara self-energy
!! functions (bare) only.
!!
  subroutine cal_sigma()
     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh

     use context, only : sigdc, sigma

     implicit none

!! local variables
     ! loop index for frequency mesh
     integer :: m

     ! loop index for spins
     integer :: s

     ! loop index for impurity sites
     integer :: t

!! [body

     ! check working axis.
     ! now only the Matsubara frequency axis is supported.
     call s_assert2(axis == 1, 'axis is wrong')

     ! loop over quantum impurities, spins, and frequency points.
     !
     ! substract the double counting terms: new sigma = sigma - sigdc.
     do t=1,nsite
         do s=1,nspin
             do m=1,nmesh
                 sigma(:,:,m,s,t) = sigma(:,:,m,s,t) - sigdc(:,:,s,t)
             enddo ! over m={1,nmesh} loop
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_sigma

!!
!! @sub cal_fermi
!!
!! try to determine the fermi level. in addition, the lattice occupancy
!! will be calculated at the same time.
!!
  subroutine cal_fermi(occup)
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nkpt, nspin
     use control, only : nmesh

     use context, only : qbnd

     implicit none

!! external arguments
     ! lattice occupancy
     real(dp), intent(out) :: occup

!! local variables
     ! desired charge density
     real(dp) :: ndens

     ! status flag
     integer  :: istat

     ! dummy array, used to save the eigenvalues of H + \Sigma(i\omega_n)
     complex(dp), allocatable :: eigs(:,:,:,:)

     ! dummy array, used to save the eigenvalues of H + \Sigma(oo)
     complex(dp), allocatable :: einf(:,:,:)

!! [body

     ! allocate memory
     allocate(eigs(qbnd,nmesh,nkpt,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_fermi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(einf(qbnd,nkpt,nspin),       stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_fermi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! calculate the nominal charge density according to the raw
     ! dft eigenvalues.
     call cal_nelect(ndens); occup = ndens

     ! construct H + \Sigma, then diagonalize it to obtain the
     ! dft + dmft eigenvalues.
     call cal_eigsys(eigs, einf)

     ! search the fermi level using bisection algorithm.
     ! the global variable `fermi` will be updated within `dichotomy()`.
     call dichotomy(ndens, eigs, einf)

     ! deallocate memory
     if ( allocated(eigs) ) deallocate(eigs)
     if ( allocated(einf) ) deallocate(einf)

!! body]

     return
  end subroutine cal_fermi

!!
!! @sub cal_eimps
!!
!! try to calculate local energy levels for all impurity sites. here,
!! eimps is defined as \sum_k \epsilon_{n,k} - \mu.
!!
  subroutine cal_eimps()
     use constants, only : dp, mystd
     use constants, only : czero

     use mmpi, only : mp_barrier
     use mmpi, only : mp_allreduce

     use control, only : nkpt, nspin
     use control, only : nsite
     use control, only : fermi
     use control, only : myid, master, nprocs

     use context, only : i_wnd
     use context, only : qdim
     use context, only : ndim
     use context, only : kwin
     use context, only : weight
     use context, only : enk
     use context, only : eimps

     implicit none

!! local variables
     ! loop index for spins
     integer :: s

     ! loop index for k-points
     integer :: k

     ! index for impurity sites
     integer :: t

     ! number of dft bands for given k-point and spin
     integer :: cbnd

     ! number of correlated orbitals for given impurity site
     integer :: cdim

     ! band window: start index and end index for bands
     integer :: bs, be

     ! status flag
     integer :: istat

     ! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

     ! dummy array, used to build site-dependent impurity level
     complex(dp), allocatable :: Xe(:,:)

     ! dummy array, used to perform mpi reduce operation for eimps
     complex(dp), allocatable :: eimps_mpi(:,:,:,:)

!! [body

     ! allocate memory
     allocate(Xe(qdim,qdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_eimps','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(eimps_mpi(qdim,qdim,nspin,nsite), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_eimps','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! reset cbnd and cdim. they will be updated later.
     ! cbnd should be k-dependent and cdim should be impurity-dependent.
     cbnd = 0
     cdim = 0

     ! reset eimps
     eimps = czero
     eimps_mpi = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate eimps for', nsite, 'sites'
         write(mystd,'(4X,a,2X,i4,2X,a)') 'add contributions from', nkpt, 'kpoints'
     endif ! back if ( myid == master ) block

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

             ! evaluate band window for the current k-point and spin.
             !
             ! i_wnd(t) returns the corresponding band window for given
             ! impurity site t. see remarks in cal_nelect().
             t = 1 ! t is fixed to 1
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

             ! determine cbnd
             cbnd = be - bs + 1

             ! provide some useful information
             write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
             write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
             write(mystd,'(2X,a,3i3)',advance='no') 'window: ', bs, be, cbnd
             write(mystd,'(2X,a,i2)') 'proc: ', myid

             ! allocate memory
             allocate(Em(cbnd),      stat = istat)
             allocate(Hm(cbnd,cbnd), stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_eimps','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

             ! evaluate `Em`, which is the eigenvalues substracted
             ! by the fermi level.
             Em = enk(bs:be,k,s) - fermi

             ! convert `Em` to diagonal matrix `Hm`
             call s_diag_z(cbnd, Em, Hm)

             ! project effective hamiltonian from the Kohn-Sham basis
             ! to the local basis, and then sum it up.
             do t=1,nsite
                 Xe = czero
                 cdim = ndim(t)
                 call one_psi_chi(cbnd, cdim, k, s, t, Hm, Xe(1:cdim,1:cdim))
                 eimps(:,:,s,t) = eimps(:,:,s,t) + Xe * weight(k)
             enddo ! over t={1,nsite} loop

             ! deallocate memory
             if ( allocated(Em) ) deallocate(Em)
             if ( allocated(Hm) ) deallocate(Hm)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! collect data from all mpi processes
# if defined (MPI)
     !
     call mp_barrier()
     !
     call mp_allreduce(eimps, eimps_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     eimps_mpi = eimps

# endif /* MPI */

     ! renormalize the impurity levels
     eimps = eimps_mpi / float(nkpt)

     ! deallocate memory
     if ( allocated(Xe) ) deallocate(Xe)
     if ( allocated(eimps_mpi) ) deallocate(eimps_mpi)

!! body]

     return
  end subroutine cal_eimps

!!
!! @sub cal_eimpx
!!
!! try to calculate local energy levels for all impurity sites. here,
!! eimpx is equal to eimps - sigdc.
!!
  subroutine cal_eimpx()
     use control, only : nspin
     use control, only : nsite

     use context, only : ndim
     use context, only : eimps, eimpx
     use context, only : sigdc

     implicit none

!! local variables
     ! index for impurity sites
     integer :: t

     ! loop index for spins
     integer :: s

     ! number of correlated orbitals for given impurity site
     integer :: cdim

!! [body

     ! substract the double counting terms from eimps to build eimpx
     do t=1,nsite
         do s=1,nspin
             cdim = ndim(t)
             eimpx(1:cdim,1:cdim,s,t) = eimps(1:cdim,1:cdim,s,t) - sigdc(1:cdim,1:cdim,s,t)
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_eimpx

!!
!! @sub cal_green
!!
!! try to calculate local green's function for all the impurity sites.
!!
  subroutine cal_green()
     use constants, only : dp
     use constants, only : czero, czi
     use constants, only : mystd

     use mmpi, only : mp_barrier
     use mmpi, only : mp_allreduce

     use control, only : nkpt, nspin
     use control, only : nsite
     use control, only : nmesh
     use control, only : myid, master, nprocs

     use context, only : i_wnd
     use context, only : qdim
     use context, only : ndim
     use context, only : kwin
     use context, only : weight
     use context, only : green

     implicit none

!! local variables
     ! loop index for spin
     integer :: s

     ! loop index for k-points
     integer :: k

     ! loop index for impurity sites
     integer :: t

     ! number of dft bands for given k-point and spin
     integer :: cbnd

     ! number of correlated orbitals for given impurity site
     integer :: cdim

     ! band window: start index and end index for bands
     integer :: bs, be

     ! status flag
     integer :: istat

     ! dummy array: for self-energy function (upfolded to Kohn-Sham basis)
     complex(dp), allocatable :: Sk(:,:,:)
     complex(dp), allocatable :: Xk(:,:,:)

     ! dummy array: for lattice green's function
     complex(dp), allocatable :: Gk(:,:,:)

     ! dummy array: for local green's function
     complex(dp), allocatable :: Gl(:,:,:)

     ! dummy array: used to perform mpi reduce operation for green
     complex(dp), allocatable :: green_mpi(:,:,:,:,:)

!! [body

     ! allocate memory for Gl and green_mpi
     allocate(Gl(qdim,qdim,nmesh), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_green','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(green_mpi(qdim,qdim,nmesh,nspin,nsite), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_green','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! reset cbnd and cdim. they will be updated later.
     ! cbnd should be k-dependent and cdim should be impurity-dependent.
     cbnd = 0
     cdim = 0

     ! reset green
     green = czero
     green_mpi = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate green for', nsite, 'sites'
         write(mystd,'(4X,a,2X,i4,2X,a)') 'add contributions from', nkpt, 'kpoints'
     endif ! back if ( myid == master ) block

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     ! loop over spins and k-points
     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

             ! evaluate band window for the current k-point and spin.
             !
             ! i_wnd(t) returns the corresponding band window for given
             ! impurity site t. see remarks in cal_nelect() for more details.
             t = 1 ! t is fixed to 1
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

             ! determine cbnd
             cbnd = be - bs + 1

             ! provide some useful information
             write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
             write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
             write(mystd,'(2X,a,3i3)',advance='no') 'window: ', bs, be, cbnd
             write(mystd,'(2X,a,i2)') 'proc: ', myid

             ! allocate memories Sk, Xk, and Gk.
             ! their sizes are k-dependent.
             allocate(Sk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Xk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Gk(cbnd,cbnd,nmesh), stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_green','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

             ! build self-energy function, and then upfold it into
             ! Kohn-Sham basis. Sk should contain contributions from
             ! all impurity sites.
             Sk = czero
             do t=1,nsite
                 Xk = czero ! reset Xk
                 cdim = ndim(t)
                 call cal_sl_sk(cdim, cbnd, k, s, t, Xk)
                 Sk = Sk + Xk
             enddo ! over t={1,nsite} loop

             ! calculate lattice green's function
             call cal_sk_gk(cbnd, bs, be, k, s, Sk, Gk)

             ! downfold the lattice green's function to obtain local
             ! green's function, then we have to perform k-summation.
             do t=1,nsite
                 Gl = czero
                 cdim = ndim(t)
                 call cal_gk_gl(cbnd, cdim, k, s, t, Gk, Gl(1:cdim,1:cdim,:))
                 green(:,:,:,s,t) = green(:,:,:,s,t) + Gl * weight(k)
             enddo ! over t={1,nsite} loop

             ! deallocate memories
             if ( allocated(Sk) ) deallocate(Sk)
             if ( allocated(Xk) ) deallocate(Xk)
             if ( allocated(Gk) ) deallocate(Gk)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! collect data from all mpi processes
# if defined (MPI)
     !
     call mp_barrier()
     !
     call mp_allreduce(green, green_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     green_mpi = green

# endif /* MPI */

     ! renormalize local green's function
     green = green_mpi / float(nkpt)

     ! deallocate memory
     if ( allocated(Gl) ) deallocate(Gl)
     if ( allocated(green_mpi) ) deallocate(green_mpi)

!! body]

     return
  end subroutine cal_green

!!
!! @sub cal_weiss
!!
!! try to calculate local weiss's function for all impurity sites.
!!
  subroutine cal_weiss()
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nspin
     use control, only : nsite, nmesh
     use control, only : myid, master

     use context, only : ndim
     use context, only : sigma
     use context, only : green
     use context, only : weiss

     implicit none

!! local variables
     ! loop index for impurity sites
     integer :: t

     ! loop index for spins
     integer :: s

     ! loop index for frequency mesh
     integer :: m

     ! number of correlated orbitals for given impurity site
     integer :: cdim

     ! status flag
     integer :: istat

     ! dummy array: for local green's function
     complex(dp), allocatable :: Gl(:,:)

!! [body

     ! reset weiss
     weiss = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate weiss for', nsite, 'sites'
     endif ! back if ( myid == master ) block

!
! remarks:
!
! try to calculate bath weiss's function using the following equation:
!
!     G^{-1}_{0} = G^{-1} + \Sigma
!
! please be aware that the double counting terms have been substracted
! from the self-energy function. see subroutine cal_sigma().
!

     ! loop over quantum impurities
     SITE_LOOP: do t=1,nsite
         ! get size of orbital space
         cdim = ndim(t)

         ! allocate memory
         allocate(Gl(cdim,cdim), stat = istat)
         !
         if ( istat /= 0 ) then
             call s_print_error('cal_weiss','can not allocate enough memory')
         endif ! back if ( istat /= 0 ) block
         !
         Gl = czero

         ! loop over spins and frequency mesh
         SPIN_LOOP: do s=1,nspin
             MESH_LOOP: do m=1,nmesh

                 ! copy local green's function to Gl
                 Gl = green(1:cdim,1:cdim,m,s,t)

                 ! inverse local green's function.
                 ! now Gl is G^{-1}.
                 call s_inv_z(cdim, Gl)

                 ! plus the self-energy function.
                 ! now Gl is G^{-1} + \Sigma.
                 Gl = Gl + sigma(1:cdim,1:cdim,m,s,t)

                 ! inverse it again to obtain bath weiss's function.
                 ! now Gl is G_0.
                 call s_inv_z(cdim, Gl)

                 ! save the final resuls to weiss
                 weiss(1:cdim,1:cdim,m,s,t) = Gl

             enddo MESH_LOOP ! over m={1,nmesh} loop
         enddo SPIN_LOOP ! over s={1,nspin} loop

         ! deallocate memory
         if ( allocated(Gl) ) deallocate(Gl)

     enddo SITE_LOOP ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_weiss

!!
!! @sub cal_delta
!!
!! try to calculate hybridization function for all impurity sites.
!!
  subroutine cal_delta()
     use constants, only : dp, mystd
     use constants, only : czi, czero

     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh
     use control, only : myid, master

     use context, only : ndim
     use context, only : fmesh
     use context, only : eimps
     use context, only : sigma
     use context, only : green
     use context, only : delta

     implicit none

!! local variables
     ! index for impurity sites
     integer :: t

     ! loop index for spin
     integer :: s

     ! loop index for frequency mesh
     integer :: m

     ! number of correlated orbitals for given impurity site
     integer :: cdim

     ! status flag
     integer :: istat

     ! dummy variables
     complex(dp) :: caux

     ! dummy arrays
     ! for identity matrix
     complex(dp), allocatable :: Im(:,:)

     ! for local impurity levels matrix
     complex(dp), allocatable :: Em(:,:)

     ! for green's function matrix
     complex(dp), allocatable :: Tm(:,:)

     ! for self-energy function matrix
     complex(dp), allocatable :: Sm(:,:)

!! [body

     ! reset delta
     delta = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate delta for', nsite, 'sites'
     endif ! back if ( myid == master ) block

     ! loop over quantum impurities
     SITE_LOOP: do t=1,nsite

         ! determine dimensional parameter
         cdim = ndim(t)

         ! allocate memory
         allocate(Im(cdim,cdim), stat = istat)
         allocate(Em(cdim,cdim), stat = istat)
         allocate(Tm(cdim,cdim), stat = istat)
         allocate(Sm(cdim,cdim), stat = istat)
         !
         if ( istat /= 0 ) then
             call s_print_error('cal_delta','can not allocate enough memory')
         endif ! back if ( istat /= 0 ) block

         ! loop over spins and frequency meshes
         SPIN_LOOP: do s=1,nspin
             MESH_LOOP: do m=1,nmesh

                 ! build identity matrix
                 call s_identity_z(cdim, Im)

                 ! get frequency point.
                 ! note that the fermi level (chemical potential) is
                 ! already included in the impurity levels `eimps`.
                 ! so here we just ignore the fermi level.
                 if ( axis == 1 ) then
                     caux = czi * fmesh(m)
                 else
                     caux = fmesh(m)
                 endif ! back if ( axis == 1 ) block

                 ! get local impurity levels.
                 ! the local impurity levels are actually equal to
                 !
                 !     \sum e_{nk} - \mu.
                 !
                 ! see cal_eimps() subroutine for more details.
                 Em = eimps(1:cdim,1:cdim,s,t)

                 ! calculate G^{-1}
                 Tm = green(1:cdim,1:cdim,m,s,t)
                 call s_inv_z(cdim, Tm)

                 ! get self-energy function.
                 ! be aware that the double counting terms have been
                 ! removed from the self-energy functions. see the
                 ! cal_sigma() subroutine for more details.
                 Sm = sigma(1:cdim,1:cdim,m,s,t)

                 ! assemble the hybridization function.
                 ! actually, Tm + Sm is G^{-1}_0.
                 ! please see cal_weiss() subroutine for more details.
                 delta(1:cdim,1:cdim,m,s,t) = caux * Im - Em - Tm - Sm

             enddo MESH_LOOP ! over m={1,nmesh} loop
         enddo SPIN_LOOP ! over s={1,nspin} loop

         ! deallocate memory
         if ( allocated(Im) ) deallocate(Im)
         if ( allocated(Em) ) deallocate(Em)
         if ( allocated(Tm) ) deallocate(Tm)
         if ( allocated(Sm) ) deallocate(Sm)

     enddo SITE_LOOP ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_delta

!!
!! @sub cal_gamma
!!
!! try to calculate correlation-induced correction for density matrix.
!!
  subroutine cal_gamma(ecorr)
     use constants, only : dp

     use control, only : nkpt, nspin

     use context, only : qbnd
     use context, only : gamma

     implicit none

!! external arguments
     ! correction to band energy
     real(dp), intent(out) :: ecorr

!! local variables
     ! status flag
     integer  :: istat

     ! dft + dmft density matrix
     complex(dp), allocatable :: kocc(:,:,:,:)

!! [body

     ! allocate memory
     allocate(kocc(qbnd,qbnd,nkpt,nspin), stat = istat)
     !
     if ( istat /= 0 ) then
         call s_print_error('cal_gamma','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! calculate new density matrix based
     ! on the dft + dmft eigenvalues.
     call cal_denmat(kocc)

     ! calculate the difference between dft
     ! and dft + dmft density matrices.
     call correction(kocc, gamma, ecorr)

     ! deallocate memory
     if ( allocated(kocc) ) deallocate(kocc)

!! body]

     return
  end subroutine cal_gamma

!!
!! @sub cal_green_tetra
!!
  subroutine cal_green_tetra()
     use constants, only : dp
     use constants, only : czero, czi
     use constants, only : mystd

     use control, only : nmesh
     use control, only : nsite
     use control, only : nkpt, nspin
     use control, only : fermi

     use context, only : qbnd
     use context, only : ndim
     use context, only : i_wnd, kwin
     use context, only : fmesh
     use context, only : green

     implicit none

!! local variables
     integer :: istat
     integer :: m
     integer :: k, s
     integer :: t
     integer :: bs, be
     integer :: cbnd
     integer :: cdim

! brillouin zone integration weight
     complex(dp), allocatable :: wtet(:,:)

! complex eigenvalues
     complex(dp), allocatable :: zenk(:,:)

! left eigenvectors, A^{L}
     complex(dp), allocatable :: zevl(:,:,:)

! right eigenvectors, A^{R}
     complex(dp), allocatable :: zevr(:,:,:)

! effective hamiltonian: hdmf = hamk + sigw(x)
     complex(dp), allocatable :: hdmf(:,:,:)

     ! dummy array: for self-energy function (upfolded to Kohn-Sham basis)
     complex(dp), allocatable :: Hk(:,:)
     complex(dp), allocatable :: Sk(:,:)
     complex(dp), allocatable :: Xk(:,:)

     ! dummy array: used to perform mpi reduce operation for green
     complex(dp), allocatable :: green_mpi(:,:,:,:,:)
     complex(dp), allocatable :: gk(:,:)
     allocate(gk(qbnd,qbnd,nkpt), stat = istat)

     allocate(wtet(qbnd,nkpt),       stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('wann_dmft_core2','can not allocate enough memory')
     endif

     allocate(zenk(qbnd,nkpt),       stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('wann_dmft_core2','can not allocate enough memory')
     endif

     allocate(zevl(qbnd,qbnd,nkpt),  stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('wann_dmft_core2','can not allocate enough memory')
     endif

     allocate(zevr(qbnd,qbnd,nkpt),  stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('wann_dmft_core2','can not allocate enough memory')
     endif

     allocate(hdmf(qbnd,qbnd,nkpt),  stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('wann_dmft_core2','can not allocate enough memory')
     endif

     ! reset green
     green = czero

     do m=1,nmesh

         do s=1,nspin
             do k=1,nkpt

                 ! evaluate band window for the current k-point and spin.
                 !
                 ! i_wnd(t) returns the corresponding band window for given
                 ! impurity site t. see remarks in cal_nelect() for more details.
                 t = 1 ! t is fixed to 1
                 bs = kwin(k,s,1,i_wnd(t))
                 be = kwin(k,s,2,i_wnd(t))

                 ! determine cbnd
                 cbnd = be - bs + 1

                 ! allocate memories Sk, Xk.
                 ! their sizes are k-dependent.
                 allocate(Hk(cbnd,cbnd), stat = istat)
                 allocate(Sk(cbnd,cbnd), stat = istat)
                 allocate(Xk(cbnd,cbnd), stat = istat)
                 !
                 if ( istat /= 0 ) then
                     call s_print_error('cal_green','can not allocate enough memory')
                 endif ! back if ( istat /= 0 ) block

                 ! build self-energy function, and then upfold it into
                 ! Kohn-Sham basis. Sk should contain contributions from
                 ! all impurity sites.
                 Sk = czero
                 do t=1,nsite
                     Xk = czero ! reset Xk
                     cdim = ndim(t)
                     call cal_sl_sk_T(cdim, cbnd, k, s, m, t, Xk)
                     Sk = Sk + Xk
                 enddo ! over t={1,nsite} loop

                 call cal_sk_hk_T(cbnd, bs, be, k, s, Sk, Hk)

                 hdmf(1:cbnd,1:cbnd,k) = Hk 

                 ! deallocate memories
                 if ( allocated(Hk) ) deallocate(Hk)
                 if ( allocated(Sk) ) deallocate(Sk)
                 if ( allocated(Xk) ) deallocate(Xk)

             enddo
         enddo

         call wann_diag_hamk3(qbnd, nkpt, hdmf, zenk, zevl, zevr)
         call wann_tetra_weight2(qbnd, czi * fmesh(m) + fermi , zenk, wtet)
         call wann_dmft_ksum2(qbnd, nkpt, wtet, zevl, zevr, gk)


     enddo

     deallocate(wtet)

     deallocate(zenk)
     deallocate(zevl)
     deallocate(zevr)

     deallocate(hdmf)

     return
  end subroutine cal_green_tetra

!>>> perform brillouin zone integration by analytical tetrahedron method
! to calculate the lattice green's function
  subroutine wann_dmft_ksum2(nwan, nkpt, wtet, zevl, zevr, gk)
     use constants, only : dp, czero

     implicit none

! external arguments
     integer, intent(in) :: nwan
     integer, intent(in) :: nkpt

! tetrahedron integration weight
     complex(dp), intent(in)  :: wtet(nwan,nkpt)

! left eigenvectors
     complex(dp), intent(in)  :: zevl(nwan,nwan,nkpt)

! right eigenvectors
     complex(dp), intent(in)  :: zevr(nwan,nwan,nkpt)

! lattice green's function
     complex(dp), intent(out) :: gk(nwan,nwan,nkpt)

! local variables
! loop index for Wannier orbitals
     integer :: iwan
     integer :: jwan
     integer :: kwan

! loop index for k-points
     integer :: ikpt

! dummy complex variables
     complex(dp) :: caux

! G_{loc}(j, i) = \sum_{ik} \sum_{k} A^{R}(j, k, ik) . W(k, ik) . A^{L}(k, i, ik)
! it is important to add up the contributions of every k-points and bands
     DMFT_WANN_LOOP1: do iwan=1,nwan               ! loop over Wannier orbitals
         DMFT_WANN_LOOP2: do jwan=1,nwan           ! loop over Wannier orbitals

             DMFT_KPNT_LOOP1: do ikpt=1,nkpt       ! loop over k-points
                 caux = czero

                 DMFT_WANN_LOOP3: do kwan=1,nwan   ! loop over Wannier orbitals
                     caux = caux + wtet(kwan,ikpt) * zevr(jwan,kwan,ikpt) * zevl(kwan,iwan,ikpt)
                 enddo DMFT_WANN_LOOP3 ! over kwan={1,n1wan} loop

                 gk(jwan,iwan,ikpt) = caux
             enddo DMFT_KPNT_LOOP1 ! over ikpt={1,nkpt} loop

         enddo DMFT_WANN_LOOP2 ! over jwan={1,n1wan} loop
     enddo DMFT_WANN_LOOP1 ! over iwan={1,n1wan} loop

     return
  end subroutine wann_dmft_ksum2

!>>> compute tetrahedron integrated weights for brillouin zone integration.
! only lambin-vigneron algorithm is implemented.
! note: it is called by wann_dmft_core1() and wann_dmft_core2() subroutines,
! used to calculate lattice green's function
  subroutine wann_tetra_weight2(nwan, z, zenk, weight)
     use constants
     use control, only : nkpt, ntet
     use context, only : tetra

     use mtetra

     implicit none

! external arguments
     integer, intent(in) :: nwan

! current complex energy
     complex(dp), intent(in)  :: z

! eigenvalues of hamiltonian
     complex(dp), intent(in)  :: zenk(nwan,nkpt)

! integration weights
     complex(dp), intent(out) :: weight(nwan,nkpt)

! local variables
! loop index for Wannier orbitals
     integer :: iwan

! loop index for tetrahedron
     integer :: itet

! loop index for k-points
     integer :: ikpt

! loop index for corner of tetrahedron
     integer :: iccc

! mass of tetrahedron
     integer, save :: mtet = 0

! energy at the corner of tetrahedron
     complex(dp) :: zc(4)

! weight at the corner of tetrahedron
     complex(dp) :: zw(4)

! build mtet for the first running
     if ( mtet == 0 ) then
         mtet = sum( tetra(:,5) )
     endif

! initialize weight
     weight = czero

     WEIGHT_TETRA_LOOP: do itet=1,ntet
         WEIGHT_ORBIT_LOOP: do iwan=1,nwan

! sets the vectors of corner energies to zc
             do iccc=1,4
                 ikpt = tetra(itet,iccc)
                 zc(iccc) = zenk(iwan,ikpt)
             enddo ! over iccc={1,4} loop

! actually calculates weights for 4 corners of one tetrahedron
             call tetra_lambin_weight(z, zc, zw)

! stores weights for irreducible k-points
             do iccc=1,4
                 ikpt = tetra(itet,iccc)
                 weight(iwan,ikpt) = weight(iwan,ikpt) + zw(iccc) * real( tetra(itet,5) )
             enddo ! over iccc={1,4} loop

         enddo WEIGHT_ORBIT_LOOP ! over iwan={1,nwan} loop
     enddo WEIGHT_TETRA_LOOP ! over itet={1,ntet} loop

! normalize properly
     weight = weight / real(mtet)

     return
  end subroutine wann_tetra_weight2

!!
!! @sub cal_sl_sk_T
!!
!! try to upfold the self-energy function from local basis to Kohn-Sham
!! basis. here, we don't care whether the double counting terms have been
!! substracted from the self-energy functions.
!!
  subroutine cal_sl_sk_T(cdim, cbnd, k, s, m, t, Sk)
     use constants, only : dp
     use constants, only : czero

     use context, only : sigma

     implicit none

!! external arguments
     ! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

     ! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

     ! index for k-points
     integer, intent(in) :: k

     ! index for spin
     integer, intent(in) :: s

     integer, intent(in) :: m

     ! index for impurity sites
     integer, intent(in) :: t

     ! self-energy function in Kohn-Sham basis
     complex(dp), intent(out) :: Sk(cbnd,cbnd)

!! local variables
     ! status flag
     integer :: istat

     ! dummy array: for local self-energy function
     complex(dp), allocatable :: Sl(:,:)

!! [body

     ! allocate memory
     allocate(Sl(cdim,cdim), stat = istat)
     !
     if ( istat /= 0 ) then
         call s_print_error('cal_sl_sk_T','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! we use Sl to store parts of sigma.
     ! note that actually sigdc has been substracted from sigma before
     ! hand. please see cal_sigma() for more details.
     Sl = sigma(1:cdim,1:cdim,m,s,t)

     ! upfolding: Sl (local basis) -> Sk (Kohn-Sham basis)
     call one_chi_psi(cdim, cbnd, k, s, t, Sl, Sk)

     ! deallocate memory
     if ( allocated(Sl) ) deallocate(Sl)

!! body]

     return
  end subroutine cal_sl_sk_T

!!
!! @sub cal_sk_hk_T
!!
!! try to build H(k) + \Sigma(i\omega_n). here, \Sigma should contain
!! contributions from all impurity sites. so, only when nsite = 1, we
!! can use the output of cal_sl_sk() as the input of this subroutine.
!!
  subroutine cal_sk_hk_T(cbnd, bs, be, k, s, Sk, Hk)
     use constants, only : dp

     use context, only : enk

     implicit none

!! external arguments
     ! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

     ! band window: start index and end index for bands
     integer, intent(in) :: bs, be

     ! index for k-points
     integer, intent(in) :: k

     ! index for spin
     integer, intent(in) :: s

     ! self-energy function at Kohn-Sham basis
     complex(dp), intent(in)  :: Sk(cbnd,cbnd)

     ! effective hamiltonian at given k-point and spin
     complex(dp), intent(out) :: Hk(cbnd,cbnd)

!! local variables
     ! status flag
     integer :: istat

     ! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

!! [body

     ! allocate memory
     allocate(Em(cbnd),      stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_hk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Hm(cbnd,cbnd), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_hk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! evaluate Em, which is just some dft eigenvalues
     Em = enk(bs:be,k,s)

     ! convert `Em` to diagonal matrix `Hm`
     call s_diag_z(cbnd, Em, Hm)

     ! combine `Hm` and `Sk` to build the effective hamiltonian
     Hk = Hm + Sk

     ! deallocate memory
     if ( allocated(Em) ) deallocate(Em)
     if ( allocated(Hm) ) deallocate(Hm)

!! body]

     return
  end subroutine cal_sk_hk_T

!>>> diagonalize general hamiltonian, return eigenvalues and eigenvectors
  subroutine wann_diag_hamk3(nwan, nkpt, hamk, eval, evecl, evecr)
     use constants, only : dp

     implicit none

! external arguments
! number of Wannier orbitals
     integer, intent(in) :: nwan

! number of k-points
     integer, intent(in) :: nkpt

! hamiltonian matrix
     complex(dp), intent(in)  :: hamk(nwan,nwan,nkpt)

! eigenvalues
     complex(dp), intent(out) :: eval(nwan,nkpt)

! left  eigenvectors, A^{L}
     complex(dp), intent(out) :: evecl(nwan,nwan,nkpt)

! right eigenvectors, A^{R}
     complex(dp), intent(out) :: evecr(nwan,nwan,nkpt)

! local variables
! loop index for k-points
     integer  :: ikpt

! dummy integer variables for lapack call
     integer  :: info
     integer  :: lwork

! working space for lapack call
     real(dp) :: rwork(2*nwan)

! working space for lapack call
     complex(dp) :: work(4*nwan)

! dummy eigenvalues and eigenvectors array
     complex(dp) :: eig(nwan)
     complex(dp) :: evl(nwan,nwan)
     complex(dp) :: evr(nwan,nwan)

! dummy hamiltonian matrix on entry, on exit, A contains intermediate values
     complex(dp) :: A(nwan,nwan)

! setup lwork
     lwork = 4 * nwan

     HAMK_DIAG_LOOP: do ikpt=1,nkpt

! copy hamk to A
         A = hamk(:,:,ikpt)

! call lapack subroutine zgeev() to diagonalize a general matrix
         call zgeev('V', 'V', nwan, A, nwan, eig, evl, nwan, evr, nwan, work, lwork, rwork, info)

! transpose conjugate left eigenvectors
         evl = transpose( dconjg( evl ) )

! make degenerate eigenvectors orthogonal
         call wann_orth_eigsys(nwan, eig, evl, evr)

! normalize eigenvectors to ensure A^{L} . A^{R} = I by Lowdin transformation
         call wann_norm_eigsys(nwan, evl, evr)

! deal with the causality of eigenvalues
         call wann_caus_eigsys(nwan, eig)

! sort final eigenvalues and corresponding eigenvectors
         call wann_sort_eigsys(nwan, eig, evl, evr)

! copy eigenvalues and eigenvectors
         eval(:,ikpt) = eig
         evecl(:,:,ikpt) = evl
         evecr(:,:,ikpt) = evr

! if error occurs
         if ( info /= 0 ) then
             call s_print_error('wann_diag_hamk3','can not diagonalize hamiltonian matrix')
         endif

     enddo HAMK_DIAG_LOOP ! over ikpt={1,nkpt} loop

     return
  end subroutine wann_diag_hamk3

!>>> make degenerate eigenvectors orthogonal
  subroutine wann_orth_eigsys(nwan, zek, evl, evr)
     use constants, only : dp, eps6

     implicit none

! external arguments
! number of Wannier orbitals
     integer, intent(in) :: nwan

! eigenvalues
     complex(dp), intent(in)    :: zek(nwan)

! left  eigenvectors, A^{L}
     complex(dp), intent(inout) :: evl(nwan,nwan)

! right eigenvectors, A^{R}
     complex(dp), intent(inout) :: evr(nwan,nwan)

! local variables
! loop index for Wannier orbitals
     integer :: pwan
     integer :: qwan

! dummy variables
     complex(dp) :: pp
     complex(dp) :: pq
     complex(dp) :: qp

     do qwan=1,nwan
         do pwan=1,qwan-1
             pq = dot_product( evl(pwan,:), evr(:,qwan) )
             pp = dot_product( evl(pwan,:), evr(:,pwan) )
             if ( abs( zek(pwan) - zek(qwan) ) < eps6 .and. abs( pq ) > eps6 ) then
                 evr(:,qwan) = evr(:,qwan) - (pq / pp) * evr(:,pwan)
             endif
         enddo ! over pwan={1,qwan-1} loop

         do pwan=1,qwan-1
             qp = dot_product( evl(qwan,:), evr(:,pwan) )
             pp = dot_product( evl(pwan,:), evr(:,pwan) )
             if ( abs( zek(pwan) - zek(qwan) ) < eps6 .and. abs( qp ) > eps6 ) then
                 evl(qwan,:) = evl(qwan,:) - (qp / pp) * evl(pwan,:)
             endif
         enddo ! over pwan={1,qwan-1} loop
     enddo ! over qwan={1,nwan} loop

     return
  end subroutine wann_orth_eigsys

!>>> normalize the eigenvectors, to ensure A^{L} . A^{R} = I
  subroutine wann_norm_eigsys(nwan, evl, evr)
     use constants, only : dp, czero

     implicit none

! external arguments
! number of Wannier orbitals
     integer, intent(in) :: nwan

! left  eigenvectors, A^{L}
     complex(dp), intent(inout) :: evl(nwan,nwan)

! right eigenvectors, A^{R}
     complex(dp), intent(inout) :: evr(nwan,nwan)

! local variables
! loop index for Wannier orbitals
     integer :: iwan
     integer :: jwan

! dummy variables
     complex(dp) :: caux
     complex(dp) :: scaler(nwan)

! apply Lowdin transformation
     do iwan=1,nwan
         caux = czero
         do jwan=1,nwan
             caux = caux + evl(iwan,jwan) * evr(jwan,iwan)
         enddo ! over jwan={1,nwan} loop
         scaler(iwan) = sqrt(caux)
     enddo ! over iwan={1,nwan} loop

     do iwan=1,nwan
         evl(iwan,:) = evl(iwan,:) / scaler(iwan)
         evr(:,iwan) = evr(:,iwan) / scaler(iwan)
     enddo ! over iwan={1,nwan} loop

     return
  end subroutine wann_norm_eigsys

!>>> to check the causality of eigenvalues. they should be negative definite
  subroutine wann_caus_eigsys(nwan, zek)
     use constants, only : dp, epss

     implicit none

! external arguments
! dimension of matrices
     integer, intent(in) :: nwan

! array of eigenvalues
     complex(dp), intent(inout) :: zek(nwan)

! local variables
! loop index
     integer  :: iwan

! real and imaginary part
     real(dp) :: re
     real(dp) :: im

! to ensure the imaginary part negative
     do iwan=1,nwan
         re =  real( zek(iwan) )
         im = aimag( zek(iwan) )

         if ( im > -epss ) then
             zek(iwan) = dcmplx(re, -epss)
         endif
     enddo ! over iwan={1,nwan} loop

     return
  end subroutine wann_caus_eigsys

!>>> this routine sorts complex(dp) eigenvalues of a matrix according 
! to its real parts with the smallest in the first slot and reorders 
! the matrices of left (row) and right (column) eigenvectors in a 
! corresponding manner.
  subroutine wann_sort_eigsys(nwan, zek, evl, evr)
     use constants, only : dp

     implicit none

! external arguments
! dimension of matrices
     integer, intent(in) :: nwan

! array of eigenvalues
     complex(dp), intent(inout) :: zek(nwan)

! matrix of left  eigenvectors (row)
     complex(dp), intent(inout) :: evl(nwan,nwan)

! matrix of right eigenvectors (column)
     complex(dp), intent(inout) :: evr(nwan,nwan)

! local parameters
! maximum value
     real(dp), parameter :: maxvalue = 1000.0_dp

! local variables
! loop index
     integer  :: p
     integer  :: q
     integer  :: idx

! minimum value
     real(dp) :: minvalue

! pointer for sorted eigenvalues
     integer  :: idxarr(nwan)

! to record which eigenvalue is sorted
     logical  :: sorted(nwan)

! save the real parts of original eigenvalues
     real(dp) :: sorton(nwan)

! sorted eigenvalues
     complex(dp) :: eval(nwan)

! sorted eigenvectors
     complex(dp) :: evec(nwan,nwan)

! initialize arrays
     idxarr = 0
     sorton = dble(zek)
     sorted = .false.

! create index array
     do p=1,nwan
         minvalue = maxvalue
         do q=1,nwan
             if ( sorted(q) == .false. .and. minvalue > sorton(q) ) then
                 minvalue = sorton(q)
                 idx = q
             endif
         enddo ! over q={1,nwan} loop
         idxarr(p) = idx
         sorted(idx) = .true.
     enddo ! over p={1,nwan} loop

! permute the eigenvalues
     do p=1,nwan
         eval(p) = zek( idxarr(p) )
     enddo ! over p={1,nwan} loop
     zek = eval

! permute the right eigenvectors
     do p=1,nwan
         evec(:,p) = evr(:,idxarr(p))
     enddo ! over p={1,nwan} loop
     evr = evec

! permute the left eigenvectors
     do p=1,nwan
         evec(p,:) = evl(idxarr(p),:)
     enddo ! over p={1,nwan} loop
     evl = evec

     return
  end subroutine wann_sort_eigsys

!>>> check the eigensystem
  subroutine wann_test_eigsys(nwan, zek, evl, evr, ham)
     use constants, only : dp, cone, epss, epst

     implicit none

! external arguments
! number of Wannier orbitals
     integer, intent(in) :: nwan

! complex eigenvalues
     complex(dp), intent(in) :: zek(nwan)

! left  eigenvectors, A^{L}
     complex(dp), intent(in) :: evl(nwan,nwan)

! right eigenvectors, A^{R}
     complex(dp), intent(in) :: evr(nwan,nwan)

! hamiltoniam
     complex(dp), intent(in) :: ham(nwan,nwan)

! local variables
! loop index for Wannier orbitals
     integer :: iwan
     integer :: jwan

! number counter
     integer :: counter

! dummy matrix
     complex(dp) :: tmat(nwan,nwan)
     complex(dp) :: zmat(nwan,nwan)

! check A^{L} . A^{R} = I
!-------------------------------------------------------------------------
     tmat = matmul(evl, evr)

     counter = 0
     do iwan=1,nwan
         if ( abs( tmat(iwan,iwan) - cone ) >= epss ) then
             counter = counter + 1
         endif
     enddo ! over iwan={1,nwan} loop

! report the error
     if ( counter > 0 ) then
         call s_print_error('wann_test_eigsys','error in evl . evr = 1')
     endif

! check A^{L} . H . A^{R} = E
!-------------------------------------------------------------------------
     tmat = matmul( matmul( evl, ham ), evr )

     counter = 0
     do iwan=1,nwan
         if ( abs( tmat(iwan,iwan) - zek(iwan) ) >= epst ) then
             counter = counter + 1
         endif
     enddo ! over iwan={1,nwan} loop

! report the error
     if ( counter > 0 ) then
         call s_print_error('wann_test_eigsys','error in evl . H . evr = zek')
     endif

! check A^{R} . E . A^{L} = H
!-------------------------------------------------------------------------
     call wann_zmat_vec(nwan, zek, zmat)
     tmat = matmul( matmul( evr, zmat ), evl )

     counter = 0
     do iwan=1,nwan
         do jwan=1,nwan
             if ( abs( tmat(iwan,jwan) - ham(iwan,jwan) ) >= epst ) then
                 counter = counter + 1
             endif
         enddo ! over jwan={1,nwan} loop
     enddo ! over iwan={1,nwan} loop

! report the error
     if ( counter > 0 ) then
         call s_print_error('wann_test_eigsys','error in evr . E . evl = H')
     endif

     return
  end subroutine wann_test_eigsys

!>>> convert a complex(dp) vector to a diagonal matrix
  subroutine wann_zmat_vec(ndim, zvec, zmat)
     use constants, only : dp, czero

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! complex(dp) vector
     complex(dp), intent(in)  :: zvec(ndim)

! object matrix
     complex(dp), intent(out) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

     zmat = czero

     do i=1,ndim
         zmat(i,i) = zvec(i)
     enddo ! over i={1,ndim} loop

     return
  end subroutine wann_zmat_vec

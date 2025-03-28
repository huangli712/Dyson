!!!-----------------------------------------------------------------------
!!! project : dyson @ jacaranda
!!! program : dmft_dump_fermi
!!!           dmft_dump_eimps
!!!           dmft_dump_eimpx
!!!           dmft_dump_eigen
!!!           dmft_dump_green
!!!           dmft_dump_weiss
!!!           dmft_dump_delta
!!!           dmft_dump_gcorr
!!! source  : dmft_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 02/23/2021 by li huang (created)
!!!           03/26/2025 by li huang (last modified)
!!! purpose : provide some subroutines to write down the calculated data.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dmft_dump_fermi
!!
!! write out calculated fermi level, lattice occupancy, and correction to
!! band energy.
!!
  subroutine dmft_dump_fermi(fermi, occup, ecorr)
     use constants, only : dp
     use constants, only : mytmp

     implicit none

!! external arguments
     ! fermi level
     real(dp), intent(in) :: fermi

     ! lattice occupancy
     real(dp), intent(in) :: occup

     ! correction to band energy
     real(dp), intent(in) :: ecorr

!! [body

     ! open data file: dmft.fermi
     open(mytmp, file='dmft.fermi', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,f16.8)') '# fermi: ', fermi
     write(mytmp,'(a9,f16.8)') '# occup: ', occup
     write(mytmp,'(a9,f16.8)') '# ecorr: ', ecorr

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_fermi

!!
!! @sub dmft_dump_eimps
!!
!! write out local impurity levels, eimps. note that the double counting
!! terms have NOT been substracted from them.
!!
  subroutine dmft_dump_eimps(eimps)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : ngrp
     use control, only : nsite

     use context, only : g_imp
     use context, only : qdim
     use context, only : ndim

     implicit none

!! external arguments
     ! local impurity levels
     complex(dp), intent(in) :: eimps(qdim,qdim,nspin,ngrp)

!! local variables
     ! loop index for spins
     integer :: s

     ! loop index for groups
     integer :: t, r

     ! loop index for correlated orbitals
     integer :: p, q

!! [body

     ! open data file: dmft.eimps
     open(mytmp, file='dmft.eimps', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

     ! write separators
     write(mytmp,*)
     write(mytmp,*)

     ! write body
     do t=1,ngrp
         !
         ! get the corresponding impurity site
         r = g_imp(t)
         !
         ! it is not a correlated group
         if ( r == 0 ) CYCLE
         !
         do s=1,nspin
             !
             ! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', r, 'spin:', s, 'dims:', ndim(t)
             !
             do q=1,ndim(t)
                 do p=1,ndim(t)
                     write(mytmp,'(2i4,2f16.8)') p, q, eimps(p,q,s,t)
                 enddo ! over p={1,ndim(t)} loop
             enddo ! over q={1,ndim(t)} loop
             !
             ! write separators
             write(mytmp,*)
             write(mytmp,*)
             !
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,ngrp} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_eimps

!!
!! @sub dmft_dump_eimpx
!!
!! write out local impurity levels, eimpx. note that the double counting
!! terms have been substracted from them.
!!
  subroutine dmft_dump_eimpx(eimpx)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : ngrp
     use control, only : nsite

     use context, only : g_imp
     use context, only : qdim
     use context, only : ndim

     implicit none

!! external arguments
     ! local impurity levels. eimpx = eimps - sigdc
     complex(dp), intent(in) :: eimpx(qdim,qdim,nspin,ngrp)

!! local variables
     ! loop index for spins
     integer :: s

     ! loop index for groups
     integer :: t, r

     ! loop index for correlated orbitals
     integer :: p, q

!! [body

     ! open data file: dmft.eimpx
     open(mytmp, file='dmft.eimpx', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

     ! write separators
     write(mytmp,*)
     write(mytmp,*)

     ! write body
     do t=1,ngrp
         !
         ! get the corresponding impurity site
         r = g_imp(t)
         !
         ! it is not a correlated group
         if ( r == 0 ) CYCLE
         !
         do s=1,nspin
             !
             ! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', r, 'spin:', s, 'dims:', ndim(t)
             !
             do q=1,ndim(t)
                 do p=1,ndim(t)
                     write(mytmp,'(2i4,2f16.8)') p, q, eimpx(p,q,s,t)
                 enddo ! over p={1,ndim(t)} loop
             enddo ! over q={1,ndim(t)} loop
             !
             ! write separators
             write(mytmp,*)
             write(mytmp,*)
             !
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,ngrp} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_eimpx

!!
!! @sub dmft_dump_eigen
!!
!! write out calculated dft + dmft eigenvalues (complex numbers).
!!
  subroutine dmft_dump_eigen(eigs)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nkpt, nspin
     use control, only : nmesh

     use context, only : xbnd
     use context, only : qwin
     use context, only : fmesh

     implicit none

!! external arguments
     ! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in) :: eigs(xbnd,nmesh,nkpt,nspin)

!! local variables
     ! number of included dft bands for given k-point and spin
     integer :: cbnd

     ! band window: start index and end index for bands
     integer :: bs, be

     ! loop index for k-points
     integer :: k

     ! loop index for spins
     integer :: s

     ! loop index for frequency grid
     integer :: m

     ! loop index for bands in band window
     integer :: q

!! [body

     ! open data file: dmft.eigen
     open(mytmp, file='dmft.eigen', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,i4)') '# nkpt : ', nkpt
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# xbnd : ', xbnd

     ! write separators
     write(mytmp,*)
     write(mytmp,*)

     ! write body
     do s=1,nspin
         do k=1,nkpt

             ! evaluate band window for the current k-point and spin.
             bs = qwin(k,s,1)
             be = qwin(k,s,2)

             ! determine cbnd
             cbnd = be - bs + 1

             ! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# kpt:', k, 'spin:', s, 'cbnd:', cbnd
             !
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 !
                 do q=1,cbnd
                     write(mytmp,'(i4,2f16.8)') q, eigs(q,m,k,s)
                 enddo ! over q={1,cbnd} loop
             enddo ! over m={1,nmesh} loop

             ! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_eigen

!!
!! @sub dmft_dump_green
!!
!! write out local green's function in matsubara frequency space.
!!
  subroutine dmft_dump_green(green)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : ngrp
     use control, only : nsite
     use control, only : nmesh

     use context, only : g_imp
     use context, only : qdim
     use context, only : ndim
     use context, only : fmesh

     implicit none

!! external arguments
     ! local green's function
     complex(dp), intent(in) :: green(qdim,qdim,nmesh,nspin,ngrp)

!! local variables
     ! loop index for spins
     integer :: s

     ! loop index for frequency grid
     integer :: m

     ! loop index for groups
     integer :: t, r

     ! loop index for correlated orbitals
     integer :: p, q

!! [body

     ! open data file: dmft.green
     open(mytmp, file='dmft.green', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

     ! write separators
     write(mytmp,*)
     write(mytmp,*)

     ! write body
     do t=1,ngrp
         !
         ! get the corresponding impurity site
         r = g_imp(t)
         !
         ! it is not a correlated group
         if ( r == 0 ) CYCLE
         !
         do s=1,nspin
             !
             ! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', r, 'spin:', s, 'dims:', ndim(t)
             !
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 !
                 do q=1,ndim(t)
                     do p=1,ndim(t)
                         write(mytmp,'(2i4,2f16.8)') p, q, green(p,q,m,s,t)
                     enddo ! over p={1,ndim(t)} loop
                 enddo ! over q={1,ndim(t)} loop
             enddo ! over m={1,nmesh} loop
             !
             ! write separators
             write(mytmp,*)
             write(mytmp,*)
             !
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,ngrp} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_green

!!
!! @sub dmft_dump_weiss
!!
!! write out local weiss's function in matsubara frequency space.
!!
  subroutine dmft_dump_weiss(weiss)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : ngrp
     use control, only : nsite
     use control, only : nmesh

     use context, only : g_imp
     use context, only : qdim
     use context, only : ndim
     use context, only : fmesh

     implicit none

!! external arguments
     ! local weiss's function
     complex(dp), intent(in) :: weiss(qdim,qdim,nmesh,nspin,ngrp)

!! local variables
     ! loop index for spins
     integer :: s

     ! loop index for frequency grid
     integer :: m

     ! loop index for groups
     integer :: t, r

     ! loop index for correlated orbitals
     integer :: p, q

!! [body

     ! open data file: dmft.weiss
     open(mytmp, file='dmft.weiss', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

     ! write separators
     write(mytmp,*)
     write(mytmp,*)

     ! write body
     do t=1,ngrp
         !
         ! get the corresponding impurity site
         r = g_imp(t)
         !
         ! it is not a correlated group
         if ( r == 0 ) CYCLE
         !
         do s=1,nspin
             !
             ! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', r, 'spin:', s, 'dims:', ndim(t)
             !
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 !
                 do q=1,ndim(t)
                     do p=1,ndim(t)
                         write(mytmp,'(2i4,2f16.8)') p, q, weiss(p,q,m,s,t)
                     enddo ! over p={1,ndim(t)} loop
                 enddo ! over q={1,ndim(t)} loop
             enddo ! over m={1,nmesh} loop
             !
             ! write separators
             write(mytmp,*)
             write(mytmp,*)
             !
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,ngrp} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_weiss

!!
!! @sub dmft_dump_delta
!!
!! write out local hybridization function in matsubara frequency space.
!!
  subroutine dmft_dump_delta(delta)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : ngrp
     use control, only : nsite
     use control, only : nmesh

     use context, only : g_imp
     use context, only : qdim
     use context, only : ndim
     use context, only : fmesh

     implicit none

!! external arguments
     ! local hybridization function
     complex(dp), intent(in) :: delta(qdim,qdim,nmesh,nspin,ngrp)

!! local variables
     ! loop index for spins
     integer :: s

     ! loop index for frequency grid
     integer :: m

     ! loop index for groups
     integer :: t, r

     ! loop index for correlated orbitals
     integer :: p, q

!! [body

     ! open data file: dmft.delta
     open(mytmp, file='dmft.delta', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

     ! write separators
     write(mytmp,*)
     write(mytmp,*)

     ! write body
     do t=1,ngrp
         !
         ! get the corresponding impurity site
         r = g_imp(t)
         !
         ! it is not a correlated group
         if ( r == 0 ) CYCLE
         !
         do s=1,nspin
             !
             ! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', r, 'spin:', s, 'dims:', ndim(t)
             !
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 !
                 do q=1,ndim(t)
                     do p=1,ndim(t)
                         write(mytmp,'(2i4,2f16.8)') p, q, delta(p,q,m,s,t)
                     enddo ! over p={1,ndim(t)} loop
                 enddo ! over q={1,ndim(t)} loop
             enddo ! over m={1,nmesh} loop
             !
             ! write separators
             write(mytmp,*)
             write(mytmp,*)
             !
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,ngrp} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_delta

!!
!! @sub dmft_dump_gcorr
!!
!! write out the dmft correction for dft density matrix.
!!
  subroutine dmft_dump_gcorr(gcorr)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nkpt, nspin

     use context, only : xbnd
     use context, only : qwin
     use context, only : kmesh

     implicit none

!! external arguments
     ! local hybridization function
     complex(dp), intent(in) :: gcorr(xbnd,xbnd,nkpt,nspin)

!! local variables
     ! number of included dft bands for given k-point and spin
     integer :: cbnd

     ! band window: start index and end index for bands
     integer :: bs, be

     ! loop index for k-points
     integer :: k

     ! loop index for spins
     integer :: s

     ! loop index for bands in band window
     integer :: p, q

!! [body

     ! open data file: dmft.gcorr
     open(mytmp, file='dmft.gcorr', form='formatted', status='unknown')

     ! write parameters
     write(mytmp,'(a9,i4)') '# nkpt : ', nkpt
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# xbnd : ', xbnd

     ! write separators
     write(mytmp,*)
     write(mytmp,*)

     ! write body
     do s=1,nspin
         do k=1,nkpt

             ! evaluate band window for the current k-point and spin.
             bs = qwin(k,s,1)
             be = qwin(k,s,2)

             ! determine cbnd
             cbnd = be - bs + 1

             ! write data for given spin and k-point
             write(mytmp,'(a,i4)') '# spin:', s
             write(mytmp,'(a,i4,2X,3f16.12)') '# kpt:', k, kmesh(k,1:3)
             write(mytmp,'(3(a,i4,2X))') '# cbnd:', cbnd, 'bs:', bs, 'be:', be
             !
             do q=1,cbnd
                 do p=1,cbnd
                     write(mytmp,'(2i4,2f16.8)') p, q, gcorr(p,q,k,s)
                 enddo ! over p={1,cbnd} loop
             enddo ! over q={1,cbnd} loop

             ! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine dmft_dump_gcorr

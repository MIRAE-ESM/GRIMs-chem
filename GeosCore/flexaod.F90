#include <define.h>
!=====================================================================
!
!     module flexaod
!
!        This module calculates the optical properties of aerosols
!        such as AOD, SSA, ASM based upon MIE algorithm.
!
!                                               written by G. Curci
!                                           modified by Seungun Lee
!
!=====================================================================

module flexaod

  use dao_mod,       only : airden, bxheight, rh
  use directory_mod, only : data_dir
  use tracer_mod,    only : stt
  use tracerid_mod,  only : idtso4, idtnh4, idtnit
  use tracerid_mod,  only : idtbcpi, idtocpi, idtbcpo, idtocpo
  use tracerid_mod,  only : idtsala, idtsalc
  use tracerid_mod,  only : idtdst1, idtdst2, idtdst3, idtdst4
  use cmn_size_mod,  only : ni => iipar, nj => jjpar, nl => llpar
  use cmn_size_mod,  only : lmaster

  implicit none

  ! I/O units
  integer,parameter :: iou=10          ! Generic use

  ! Parameters
  integer,parameter :: nwl = 61        ! number input wavelengths
  integer,parameter :: nrh = 5         ! number of input RH bins
  integer,parameter :: nsinyuk = 109
  integer,parameter :: nbnd = 30
  integer,parameter :: nspecs = 6      ! Number of GEOS-Chem species
  integer,parameter :: ndust = 7       ! Number of GEOS-Chem DUST species
  integer,parameter :: npcoef = 100    ! Number of Legendre expansion coeff.
  integer,parameter :: nsizspc = nspecs + ndust - 1
  character(len=4),dimension(nspecs),parameter :: &
    species = (/ 'sulf', 'oc', 'bc', 'ssa', 'ssc', 'dust' /)
  real,dimension(6), parameter :: &
    hgf_bc_chin = (/ 1.0, 1.0, 1.0, 1.2, 1.4, 1.5 /)

  ! Module variables
  real,dimension(nbnd) :: wl_um
  integer,dimension(:),allocatable    :: rhbins
  real,dimension(:),allocatable       :: wlbins
  real,dimension(:,:),allocatable     :: q_dext, c_dext, c_dsca
  real,dimension(:,:),allocatable     :: dssalb, dasym, v_dave, r_dust
  real,dimension(:,:,:),allocatable   :: mr, mi
  real,dimension(:,:,:),allocatable   :: q_ext, r_eff, v_ave
  real,dimension(:,:,:),allocatable   :: c_ext, c_sca
  real,dimension(:,:,:),allocatable   :: ssalb, asym
  real,dimension(:,:,:),allocatable   :: specmr, specmi
  real,dimension(:,:,:,:),allocatable :: outod
  real,dimension(:,:,:,:),allocatable :: outssa, outg
  real,dimension(nsizspc)             :: rho, sigma, gam_a, gam_b
  real,dimension(:,:),allocatable     :: radius, rmin, rmax

  ! Water refractive index
  real :: mwatr, mwati

  ! misc
  integer :: ios, i, j, l, n
  integer :: ispec, irh, ibnd
  real    :: weight, dry_frac
  character(len=300) :: filename

  !====================================================================
  ! FLEXAOD begins here!
  !====================================================================

  contains

!======================================================================

  subroutine read_mie

  ! Wavelenght in microns (um)
  wl_um(1)  = (  3.846  +  3.077  ) / 2d0
  wl_um(2)  = (  3.077  +  2.500  ) / 2d0
  wl_um(3)  = (  2.500  +  2.150  ) / 2d0
  wl_um(4)  = (  2.150  +  1.942  ) / 2d0
  wl_um(5)  = (  1.942  +  1.626  ) / 2d0
  wl_um(6)  = (  1.626  +  1.299  ) / 2d0
  wl_um(7)  = (  1.299  +  1.242  ) / 2d0
  wl_um(8)  = (  1.242  +  0.7782 ) / 2d0
  wl_um(9)  = (  0.7782 +  0.6250 ) / 2d0
  wl_um(10) = (  0.6250 +  0.4415 ) / 2d0
  wl_um(11) = (  0.4415 +  0.3448 ) / 2d0
  wl_um(12) = (  0.3448 +  0.2632 ) / 2d0
  wl_um(13) = (  0.2632 +  0.2000 ) / 2d0
  wl_um(14) = ( 12.195  +  3.846  ) / 2d0
  wl_um(15) = ( 1000.0  + 28.571  ) / 2d0
  wl_um(16) = ( 28.571  + 20.000  ) / 2d0
  wl_um(17) = ( 20.000  + 15.873  ) / 2d0
  wl_um(18) = ( 15.873  + 14.286  ) / 2d0
  wl_um(19) = ( 14.286  + 12.195  ) / 2d0
  wl_um(20) = ( 12.195  + 10.204  ) / 2d0
  wl_um(21) = ( 10.204  +  9.259  ) / 2d0
  wl_um(22) = (  9.259  +  8.475  ) / 2d0
  wl_um(23) = (  8.475  +  7.194  ) / 2d0
  wl_um(24) = (  7.194  +  6.757  ) / 2d0
  wl_um(25) = (  6.757  +  5.556  ) / 2d0
  wl_um(26) = (  5.556  +  4.808  ) / 2d0
  wl_um(27) = (  4.808  +  4.444  ) / 2d0
  wl_um(28) = (  4.444  +  4.202  ) / 2d0
  wl_um(29) = (  4.202  +  3.846  ) / 2d0
  wl_um(30) = (  3.846  +  3.077  ) / 2d0

  ! Cap
  wl_um = max(0.250,wl_um)
  wl_um = min(40.00,wl_um)

  ! Allocate arrays
  call alloc_other

  !--------------------------------------------------------------------
  ! Read Optical properties
  !--------------------------------------------------------------------

  ! Message
  if (lmaster) then
  write(6,'("  o Reading optical properties ... ")')
  endif

  ! Allocate arrays
  call alloc_optics

  mr = 0.0
  mi = 0.0

      ! Hardwire RH bins
      if(nrh>0) rhbins(1) = 0
      if(nrh>1) rhbins(2) = 50
      if(nrh>2) rhbins(3) = 70
      if(nrh>3) rhbins(4) = 80
      if(nrh>4) rhbins(5) = 90
      if(nrh>5) rhbins(6) = 95
      if(nrh>6) rhbins(7) = 98
      if(nrh>7) rhbins(8) = 99

      ! Read spectral complex refractive indices
      specmr = 0.0
      specmi = 0.0
      do ispec = 1,nspecs
        do irh = 1,nrh

          if (lmaster) then
          ! Read data
          write(filename,'("opac_",A,I2.2,".dat")') &
            trim(species(ispec)),rhbins(irh)
          open(iou,file=trim(data_dir)//'/flexaod/'//filename,iostat=ios)
          if (ios /= 0) stop 'flexaod : cannot open OPAC ref. index'
          do i = 1,nwl
            read(iou,*) wlbins(i), specmr(i,irh,ispec), specmi(i,irh,ispec)
          enddo
          close(iou)
          endif

          ! Exit immediately if BC and dust (dry data only)
          if (species(ispec)=='bc' .or. species(ispec)=='dust') exit

        enddo
      enddo

#ifdef MP
#ifndef RMP
      call mpbcastr(wlbins,nwl)
      call mpbcastr(specmr,nwl*nrh*nspecs)
      call mpbcastr(specmi,nwl*nrh*nspecs)
#else
      call rmpbcastr(wlbins,nwl)
      call rmpbcastr(specmr,nwl*nrh*nspecs)
      call rmpbcastr(specmi,nwl*nrh*nspecs)
#endif
#endif

      ! Interpolation of refractive indices onto output wavelength
      do ibnd = 1, nbnd
        do ispec = 1, nspecs
          do irh = 1, nrh
            mr(irh,ispec,ibnd) = linterp( wlbins, specmr(:,irh,ispec), wl_um(ibnd) )
            mi(irh,ispec,ibnd) = linterp( wlbins, specmi(:,irh,ispec), wl_um(ibnd) )
          enddo
        enddo
      enddo

    ! Reallocate buffer arrays
    if (allocated(wlbins)) deallocate(wlbins)
    if (allocated(specmr)) deallocate(specmr)
    if (allocated(specmi)) deallocate(specmi)
    allocate(wlbins(nsinyuk))
    allocate(specmr(nsinyuk,1,1))
    allocate(specmi(nsinyuk,1,1))

    if (lmaster) then
    ! Read data
    filename = "refrac.sinyukdust"
    open(iou,file=trim(data_dir)//'/flexaod/'//filename,iostat=ios)
    if (ios /= 0) stop 'flexaod : cannot open SINYUK dust ref. index'
    do i = 1,nsinyuk
      read(iou,'(6X,F6.3,4X,F5.3,4X,E9.3)') wlbins(i), specmr(i,1,1), specmi(i,1,1)
    enddo
    close(iou)
    endif

#ifdef MP
#ifndef RMP
    call mpbcastr(wlbins,nsinyuk)
    call mpbcastr(specmr,nsinyuk)
    call mpbcastr(specmi,nsinyuk)
#else
    call rmpbcastr(wlbins,nsinyuk)
    call rmpbcastr(specmr,nsinyuk)
    call rmpbcastr(specmi,nsinyuk)
#endif
#endif

    do ibnd = 1, nbnd
      ! Interpolation of refractive indices onto output wavelength
      mr(1,nspecs,ibnd) = linterp( wlbins, specmr(:,1,1), max(wl_um(ibnd),0.3) )
      ! Change sign of imaginary index for consistency with Mie calculations
      mi(1,nspecs,ibnd) = - linterp( wlbins, specmi(:,1,1), max(wl_um(ibnd),0.3) )
    enddo

  !--------------------------------------------------------------------
  ! Read Water optical properties
  !--------------------------------------------------------------------

      mwatr = 1.33
      mwati = 0.0

  !--------------------------------------------------------------------
  ! Read size distributions
  !--------------------------------------------------------------------

  ! Message
  if (lmaster) then
  write(6,'("  o Reading size distributions ... ")')
  endif

  ! Allocate arrays
  call alloc_dist

  if (lmaster) then
  ! Read size distributions
  open(iou,file=trim(data_dir)//'/flexaod/'//"size_dist.dat",iostat=ios)
  if (ios /= 0) stop 'flexaod : cannot open size distributions'
  do i = 1,nsizspc
    read(iou,*) rho(i), radius(1,i), sigma(i), &
                gam_a(i), gam_b(i), rmin(1,i), rmax(1,i)
  enddo
  close(iou)
  endif

#ifdef MP
#ifndef RMP
  call mpbcastr(rho,nsizspc)
  call mpbcastr(radius,nrh*nsizspc)
  call mpbcastr(sigma,nsizspc)
  call mpbcastr(gam_a,nsizspc)
  call mpbcastr(gam_b,nsizspc)
  call mpbcastr(rmin,nrh*nsizspc)
  call mpbcastr(rmax,nrh*nsizspc)
#else
  call rmpbcastr(rho,nsizspc)
  call rmpbcastr(radius,nrh*nsizspc)
  call rmpbcastr(sigma,nsizspc)
  call rmpbcastr(gam_a,nsizspc)
  call rmpbcastr(gam_b,nsizspc)
  call rmpbcastr(rmin,nrh*nsizspc)
  call rmpbcastr(rmax,nrh*nsizspc)
#endif
#endif

  !====================================================================
  ! Hygroscopic growth

      if (lmaster) then
      ! Read spectral complex refractive indices
      do ispec = 1,nspecs-1

        ! Open data file
        write(filename,'("opac_",A,"_dist.dat")') trim(species(ispec))
        open(iou,file=trim(data_dir)//'/flexaod/'//filename,iostat=ios)
        if (ios /= 0) stop 'flexaod : cannot open OPAC size distributions'
        do irh = 1,nrh

          ! Read (and eventually overwrite) size dist. parameters
          read(iou,*) rmin(irh,ispec)
          read(iou,*) rmax(irh,ispec)

          ! Do not overwrite sigma if "STD" case
            read(iou,*)

          ! If BC exit now (dry data only)
          if (species(ispec)=='bc') then
            ! Do not overwrite radius if "STD" case
              read(iou,*)
            exit

          ! Other species ...
          else
            read(iou,*) radius(irh,ispec)
          endif

        enddo  ! rh bins
        close(iou)

      enddo  ! species
      endif

#ifdef MP
#ifndef RMP
  call mpbcastr(radius,nrh*nsizspc)
  call mpbcastr(rmin,nrh*nsizspc)
  call mpbcastr(rmax,nrh*nsizspc)
#else
  call rmpbcastr(radius,nrh*nsizspc)
  call rmpbcastr(rmin,nrh*nsizspc)
  call rmpbcastr(rmax,nrh*nsizspc)
#endif
#endif

      ! Hardwire BC hygroscopic growth following Chin et al., 2002
      ispec = 3
      do irh = 2,nrh
        radius(irh,ispec) = radius(1,ispec) * hgf_bc_chin(irh)
        rmin(irh,ispec) = rmin(1,ispec) * hgf_bc_chin(irh)
        rmax(irh,ispec) = rmax(1,ispec) * hgf_bc_chin(irh)
      enddo

        ! Divide Rmax by 40 for sulfate, OC (Chin et al., 2002)
        rmax(:,1) = rmax(:,1) * 0.025
        rmax(:,2) = rmax(:,2) * 0.025

        ! Multiply OC radius by 3 (Drury et al., 2010)
        radius(:,2) = radius(:,2) * 3.0

        ! Divide SSa radius by 2.45 (Jaegle et al., 2011)
        radius(:,4) = radius(:,4) / 2.45
        rmin(:,4) = 5.e-3
        rmax(:,4) = 20.

        ! Divide SSc radius by 4.36 (Jaegle et al., 2011)
        radius(:,5) = radius(:,5) / 4.36
        rmin(:,5) = 5.e-3
        rmax(:,5) = 60.

  !====================================================================

  ! RH influence on BC refractive index: mix with that of water
  ispec = 3
  do irh = 2, nrh
    dry_frac = radius(1,ispec)**3 / radius(irh,ispec)**3
    mr(irh,ispec,:) = mr(1,ispec,:)*dry_frac + mwatr*(1.-dry_frac)
    mi(irh,ispec,:) = mi(1,ispec,:)*dry_frac + mwati*(1.-dry_frac)
  enddo

  !--------------------------------------------------------------------
  ! Calculate Mie tabulated parameters
  !--------------------------------------------------------------------

    if (lmaster) then
    filename = trim(data_dir)//'/flexaod/'//'mie_tables_rrtmg'
    open(iou, file=filename, form='unformatted', status='old', iostat=ios)
    endif

#ifdef MP
#ifndef RMP
    call mpbcasti(ios,1)
#else
    call rmpbcasti(ios,1)
#endif
#endif

    if (ios .ne. 0) then
      if (lmaster) then
      write(6,'("  o Reading Mie tables failed !!! ")')
      endif
      go to 100
    endif

    if (lmaster) then
    read(iou) q_ext, c_ext, c_sca, r_eff, v_ave, ssalb, asym, &
              q_dext, c_dext, c_dsca, r_dust, v_dave, dssalb, dasym
    close(iou)
    endif

#ifdef MP
#ifndef RMP
    call mpbcastr(q_ext,nrh*(nspecs-1)*nbnd)
    call mpbcastr(c_ext,nrh*(nspecs-1)*nbnd)
    call mpbcastr(c_sca,nrh*(nspecs-1)*nbnd)
    call mpbcastr(r_eff,nrh*(nspecs-1)*nbnd)
    call mpbcastr(v_ave,nrh*(nspecs-1)*nbnd)
    call mpbcastr(ssalb,nrh*(nspecs-1)*nbnd)
    call mpbcastr(asym,nrh*(nspecs-1)*nbnd)
    call mpbcastr(q_dext,ndust*nbnd)
    call mpbcastr(c_dext,ndust*nbnd)
    call mpbcastr(c_dsca,ndust*nbnd)
    call mpbcastr(r_dust,ndust*nbnd)
    call mpbcastr(v_dave,ndust*nbnd)
    call mpbcastr(dssalb,ndust*nbnd)
    call mpbcastr(dasym,ndust*nbnd)
#else
    call rmpbcastr(q_ext,nrh*(nspecs-1)*nbnd)
    call rmpbcastr(c_ext,nrh*(nspecs-1)*nbnd)
    call rmpbcastr(c_sca,nrh*(nspecs-1)*nbnd)
    call rmpbcastr(r_eff,nrh*(nspecs-1)*nbnd)
    call rmpbcastr(v_ave,nrh*(nspecs-1)*nbnd)
    call rmpbcastr(ssalb,nrh*(nspecs-1)*nbnd)
    call rmpbcastr(asym,nrh*(nspecs-1)*nbnd)
    call rmpbcastr(q_dext,ndust*nbnd)
    call rmpbcastr(c_dext,ndust*nbnd)
    call rmpbcastr(c_dsca,ndust*nbnd)
    call rmpbcastr(r_dust,ndust*nbnd)
    call rmpbcastr(v_dave,ndust*nbnd)
    call rmpbcastr(dssalb,ndust*nbnd)
    call rmpbcastr(dasym,ndust*nbnd)
#endif
#endif

    if (lmaster) then
    write(6,'("  o Reading Mie tables done !!! ")')
    endif

    return

100 continue

    ! Message
    if (lmaster) then
    write(6,'("  o Calculating Mie tables ... ")')
    endif

    ! Mie calculations
    call mie_tab

    ! Message
    if (lmaster) then
    write(6,'("  o Mie tables done !!! ")')
    endif

    if (lmaster) then
    filename = trim(data_dir)//'/flexaod/'//'mie_tables_rrtmg'
    open(iou, file=filename, form='unformatted', status='new')
    write(iou) q_ext, c_ext, c_sca, r_eff, v_ave, ssalb, asym, &
               q_dext, c_dext, c_dsca, r_dust, v_dave, dssalb, dasym
    close(iou)
    endif

    ! Message
    if (lmaster) then
    write(6,'("  o Writing Mie tables done !!! ")')
    endif

  end subroutine read_mie

!======================================================================
! LINTERP: linear interpolation
!======================================================================
real function linterp( bins, values, outpoint )

  ! args
  real,dimension(:) :: bins, values
  real :: outpoint

  ! local
  integer :: nbins, i
  real :: weight

  ! Check we are inside the input range
  if ( outpoint<minval(bins) .or. outpoint>maxval(bins) ) then
    if (lmaster) then
    write(6,'("!!! Output point must be in the range ",F6.3," - ", &
                   F6.3," !!!")') &
           minval(bins), maxval(bins)
    endif
    stop
  endif

  ! Get number of bins
  nbins = size( bins, 1 )

  ! Test simple case: we have exact point in input
  do i = 1, nbins
    if (bins(i)==outpoint) exit
  enddo
  ! ... if YES : just copy arrays
  if (i <= nbins) then
    linterp = values(i)

  ! ... if NO : linear interpolation
  else
    do i = 2,nbins
      if (outpoint<=bins(i)) exit
    enddo
    weight = (bins(i)-outpoint) / (bins(i)-bins(i-1))
    linterp = values(i-1)*    weight + &
              values(i  )*(1.-weight)
  endif

end function linterp

!======================================================================
! MIE_TAB: Calculate RH-dependent Q and Reff at fixed RH bins
!======================================================================
subroutine mie_tab

  integer :: ispec, irh, idst
  real :: rwet, mrwet, miwet
  real*8 :: aa,bb,aa1,aa2,bb1,bb2,gam,lam,mrr,mri
  real*8 :: r1,r2,ddelt,reff,cext,csca,cbac,area,vol
  real*8 :: alb
  real*8, dimension(npcoef) :: al1
  integer :: ndistr,nk,n,np
  integer :: ip

  do ibnd = 1, nbnd

  ! Message
  if (lmaster) then
  write(6,'("   > For Wavelenth = ",F6.4," um ")') wl_um(ibnd)
  endif

  ! Non-dust species
  do ispec = 1, nspecs-1

    ! Message
    if (lmaster) then
    write(6,'("    > Tabulating species ",A," ... ")') species(ispec)
    endif

    do irh = 1, nrh

      ! Wet radius (um)
      rwet = radius(irh,ispec)

      ! Wet refractive index
      mrwet = mr(irh,ispec,ibnd)
      miwet = mi(irh,ispec,ibnd)

      ! Mie parameters
      AA=dble(rwet)
      BB=DLOG(dble(sigma(ispec)))*DLOG(dble(sigma(ispec)))
      AA1=0.0d0
      AA2=0.0d0
      BB1=0.0d0
      BB2=0.0d0
      GAM=0.0d0
      LAM=dble(wl_um(ibnd))
      MRR=dble(mrwet)
      MRI=-dble(miwet)
      NDISTR=2
      NK=100
      N=10
      NP=4
      R1=dble(rmin(irh,ispec))
      R2=dble(rmax(irh,ispec))
      DDELT=1D-5

      ! Mie calculations
      !  Cext = extinction cross-section (um^2)
      !  Csca = scattering cross-section (um^2)
      !  Cbac = backscattering cross-section (um^2 srad^-1)
      !  area = <G> = average geometric particle cross-section (um^2)
      !  vol  = <V> = average geometric particle volume (um^3)
      !  Reff = effective radius (um)
      !  alb  = Single-scattering albedo (unitless)
      !  al1  = Legendre coeffs of the phase function (unitless)
      call mie_spher(aa,bb,aa1,aa2,bb1,bb2,gam,&
                     lam,mrr,mri,ndistr,       &
                     nk,n,np,r1,r2,ddelt,      &
                     cext,csca,cbac,area,vol,  &
                     reff,alb,al1)

      ! Store Mie Table
      !  Qext = extinction efficiency (unitless)
      !  asym = asymmetry parameter (unitless)
      q_ext(irh,ispec,ibnd) = real(cext/area)
      c_ext(irh,ispec,ibnd) = real(cext)
      c_sca(irh,ispec,ibnd) = real(csca)
!      c_bac(irh,ispec,ibnd) = real(cbac)
      v_ave(irh,ispec,ibnd) = real(vol)
      ssalb(irh,ispec,ibnd) = real(alb)
      asym(irh,ispec,ibnd) = real(al1(2))/3.
      r_eff(irh,ispec,ibnd) = real(reff)

      ! debug
      !write(*,'(I1,2X,F6.4,2X,F5.3,2x,F6.4,1X,F5.3)')irh,q_ext(irh,ispec),r_eff(irh,ispec),ssalb(irh,ispec),al1(2)
    enddo
  enddo

  ! Dust species
  ! Message
  if (lmaster) then
  write(6,'("    > Tabulating species ",A," ... ")') species(6)
  endif
  do ispec = 1, ndust

    ! Dust index
    idst = nspecs - 1 + ispec

    ! Mie parameters
    AA=dble(gam_a(idst))
    BB=dble(gam_b(idst))
    AA1=0.0d0
    AA2=0.0d0
    BB1=0.0d0
    BB2=0.0d0
    GAM=0.0d0
    LAM=dble(wl_um(ibnd))
    MRR=dble(mr(1,nspecs,ibnd))
    MRI=-dble(mi(1,nspecs,ibnd))
    NDISTR=4
    NK=100
    N=10
    NP=4
    R1=dble(rmin(1,idst))
    R2=dble(rmax(1,idst))
    DDELT=1d-5

    ! Mie calculations
    !  Cext = extinction cross-section (um^2)
    !  Csca = scattering cross-section (um^2)
    !  Cbac = backscattering cross-section (um^2 srad^-1)
    !  area = <G> = average geometric particle cross-section (um^2)
    !  vol  = <V> = average geometric particle volume (um^3)
    !  Reff = effective radius (um)
    !  alb  = Single-scattering albedo (unitless)
    !  al1  = Legendre coeffs of the phase function (unitless)
    call mie_spher(aa,bb,aa1,aa2,bb1,bb2,gam,&
                   lam,mrr,mri,ndistr,       &
                   nk,n,np,r1,r2,ddelt,      &
                   cext,csca,cbac,area,vol,  &
                   reff,alb,al1)

    ! Store Mie Table
    !  Qext = extinction efficiency (unitless)
    !  asym = asymmetry parameter (unitless)
    q_dext(ispec,ibnd) = real(cext/area)
    c_dext(ispec,ibnd) = real(cext)
    c_dsca(ispec,ibnd) = real(csca)
!    c_dbac(ispec,ibnd) = real(cbac)
    v_dave(ispec,ibnd) = real(vol)
    dssalb(ispec,ibnd) = real(alb)
    dasym(ispec,ibnd) = real(al1(2))/3.
    r_dust(ispec,ibnd) = real(reff)

    ! debug
    !write(*,'(I1,2X,F6.4,2X,F5.3,2x,F5.3,8(1X,F5.3))')irh,q_dext(ispec),r_dust(ispec),dssalb(ispec),al1(1:8)
  enddo

  enddo

end subroutine mie_tab

!======================================================================
! INTERP_AOD: Calculate RH-dependent AOD from Mie table
!======================================================================
subroutine calc_aod

  integer :: i, j, l, ispec, irh, iphob, idst, r, ip
  integer :: n, icount
  real :: dz, fwet, weight
  real :: reff, qext, optdep, g, ss
  real :: scaleq, scaler, scaleod
  real :: x, y
  real :: kappa_ext, kappa_sca, totssa, totg
  real :: k_ext, k_sca, numconc
  real,dimension(nsizspc+2) :: conc_gcm3
  real,dimension(nrh) :: rw, qw, gw, ssw

  ! Initialize
  outod = 0.0
  outssa = 0.0
  outg = 0.0

!$omp parallel do collapse( 3 ) &
!$omp private( ibnd, i, j, l, ispec, irh, iphob, idst, r ) &
!$omp private( dz, fwet, weight ) &
!$omp private( reff, qext, optdep, g, ss ) &
!$omp private( scaleq, scaler, scaleod ) &
!$omp private( kappa_ext, kappa_sca, totssa, totg ) &
!$omp private( k_ext, k_sca, numconc ) &
!$omp private( conc_gcm3 ) &
!$omp private( rw, qw, gw, ssw )
  do l = 1,nl
  do j = 1,nj
  do i = 1,ni

      conc_gcm3 = 0d0
      ! Tracer concentrations (kg/kg)
      if (idtso4 .ne. 0) then
        if (idtnh4 .eq. 0) then
          conc_gcm3(1) = stt(i,j,l,idtso4) * ( 96.0 + 36.0 ) / 28.97
        else
          conc_gcm3(1) = stt(i,j,l,idtso4) * 96.0 / 28.97
        endif
      endif
      if (idtnh4 .ne. 0) then
        conc_gcm3(1) = conc_gcm3(1) + stt(i,j,l,idtnh4) * 18.0 / 28.97
      endif
      if (idtnit .ne. 0) then
        conc_gcm3(1) = conc_gcm3(1) + stt(i,j,l,idtnit) * 62.0 / 28.97
      endif
      if (idtocpi .ne. 0) then
        conc_gcm3(2) = stt(i,j,l,idtocpi) * 12.0 * 2.1 / 28.97
      endif
      if (idtbcpi .ne. 0) then
        conc_gcm3(3) = stt(i,j,l,idtbcpi) * 12.0 / 28.97
      endif
      if (idtsala .ne. 0) then
        conc_gcm3(4) = stt(i,j,l,idtsala) * 36.0 / 28.97
      endif
      if (idtsalc .ne. 0) then
        conc_gcm3(5) = stt(i,j,l,idtsalc) * 36.0 / 28.97
      endif
      if (idtocpo .ne. 0) then
        conc_gcm3(6) = stt(i,j,l,idtocpo) * 12.0 * 2.1 / 28.97
      endif
      if (idtbcpo .ne. 0) then
        conc_gcm3(7) = stt(i,j,l,idtbcpo) * 12.0 / 28.97
      endif
      if (idtdst1 .ne. 0) then
        conc_gcm3(8) = stt(i,j,l,idtdst1)*0.25 * 29.0 / 28.97
        conc_gcm3(9) = stt(i,j,l,idtdst1)*0.25 * 29.0 / 28.97
        conc_gcm3(10) = stt(i,j,l,idtdst1)*0.25 * 29.0 / 28.97
        conc_gcm3(11) = stt(i,j,l,idtdst1)*0.25 * 29.0 / 28.97
      endif
      if (idtdst2 .ne. 0) then
        conc_gcm3(12) = stt(i,j,l,idtdst2) * 29.0 / 28.97
      endif
      if (idtdst3 .ne. 0) then
        conc_gcm3(13) = stt(i,j,l,idtdst3) * 29.0 / 28.97
      endif
      if (idtdst4 .ne. 0) then
        conc_gcm3(14) = stt(i,j,l,idtdst4) * 29.0 / 28.97
      endif

      ! Aerosol concentration (g/cm3)
      conc_gcm3 = conc_gcm3 * airden(l,i,j) * 1e-3

      ! Box Height (m --> cm)
      dz = bxheight(i,j,l) * 1e2

      ! Select RH bin
      do irh = 2,nrh
        if (rh(i,j,l)<rhbins(irh)) exit
      enddo
      irh = irh - 1

      if (irh<nrh) then
        ! Interpolation weight
        weight = (rh(i,j,l)-rhbins(irh)) / (rhbins(irh+1)-rhbins(irh))
        if ( weight > 1.0d0 ) weight = 1.0d0
      endif


      ! Loop over bands
      do ibnd = 1,nbnd

        ! Reset total variables
        kappa_ext = 0.0
        kappa_sca = 0.0
        totssa = 0.0
        totg = 0.0


      ! Loop over non-dust species
      do ispec = 1,nspecs-1

        !---------------------
        ! Hygroscopic growth
        !---------------------

            ! Define wet optical quantities at OPAC RH bins
            do r = 1, nrh

              ! Wet radius
              rw(r) = r_eff(r,ispec,ibnd)

              ! Wet frac of aerosol
              fwet  = (rw(r)**3 - rw(1)**3) / rw(r)**3

              ! Wet extinction efficiency
              !qw(r) = q_ext(r,ispec)*fwet + q_ext(1,ispec)*(1.d0-fwet)
              qw(r) = q_ext(r,ispec,ibnd)

              ! Wet single-scattering albed
              !ssw(r) = ssalb(r,ispec)*fwet + ssalb(1,ispec)*(1.d0-fwet)
              ssw(r) = ssalb(r,ispec,ibnd)

              ! Wet asymmetry factor
              !gw(r) = asym(r,ispec)*fwet + asym(1,ispec)*(1.d0-fwet)
              gw(r) = asym(r,ispec,ibnd)

            enddo  ! RH bins

            ! Interpolate optical parameters
            if (irh<nrh) then

              ! Interpolate radius and Q scaling factor
              scaleq = (weight*qw(irh+1) + (1.d0-weight)*qw(irh)) / qw(1)
              reff   =  weight*rw(irh+1) + (1.d0-weight)*rw(irh)

              ! Interpolate single-scattering albedo
              ss =  weight*ssw(irh+1) + (1.d0-weight)*ssw(irh)

              ! Interpolate asymmetry parameter
              g =  weight*gw(irh+1) + (1.d0-weight)*gw(irh)

            ! Last RH bin: do not interpolate
            else

              ! Radius and Q scaling factor
              scaleq = qw(nrh) / qw(1)
              reff = rw(nrh)

              ! Single-scattering albedo
              ss =  ssw(irh)

              ! Asymmetry parameter
              g =  gw(irh)

            endif   ! RH bins

            ! Scaling parameters
            scaler = reff / rw(1)
            scaleod = scaleq * scaler * scaler
            !print *, ispec,i,j,l,rh,irh,scaler,scaleq,reff

            ! Compute box optical depth
            ! Convert Reff um --> cm
            qext = q_ext(1,ispec,ibnd)
            reff = r_eff(1,ispec,ibnd) * 1e-4
            optdep = 0.75 * qext * dz*conc_gcm3(ispec) / &
                     ( reff * rho(ispec) ) * scaleod

            ! Add hydrophobic contribution
            select case (species(ispec))
              case ('oc')
                iphob = 6
              case ('bc')
                iphob = 7
              case default
                iphob = 0
            end select
            if (iphob>0) then
              optdep = optdep + &
                       0.75 * qext * dz*conc_gcm3(iphob) / &
                       ( reff * rho(ispec) )
            endif

            ! Number concentration [#/cm3]
            numconc = conc_gcm3(ispec) / &
              (rho(ispec) * v_ave(1,ispec,ibnd)*1e-12)
            ! Extinction coefficient [km-1]
            !  Cext = extinction cross-section (um^2)
            !  Ncon = number concentration (cm^-3)
            !  Kext = extinction coefficient (km^-1)
            k_ext = c_ext(1,ispec,ibnd) * &
              numconc * scaleod*1e-3
            kappa_ext = kappa_ext + k_ext
            ! Scattering coefficient [km-1]
            k_sca = c_sca(1,ispec,ibnd) * &
              numconc * scaleod*1e-3
            kappa_sca = kappa_sca + k_sca
            ! Single Scattering Albedo
            totssa = totssa + k_ext * ss
            ! Asymmetry factor
            totg = totg + k_sca * g

            ! Hydrophobic species
            if (iphob>0) then
              ! Number concentration [#/cm3]
              numconc = conc_gcm3(iphob) / &
                (rho(ispec) * v_ave(1,ispec,ibnd)*1e-12)
              ! Extinction [km-1]
              k_ext = c_ext(1,ispec,ibnd) * &
                numconc * 1e-3
              kappa_ext = kappa_ext + k_ext
              ! Scattering [km-1]
              k_sca = c_sca(1,ispec,ibnd) * &
                numconc * 1e-3
              kappa_sca = kappa_sca + k_sca
              ! Single Scattering Albedo
              totssa = totssa + k_ext * ssalb(1,ispec,ibnd)
              ! Asymmetry factor
              totg = totg + k_sca * asym(1,ispec,ibnd)
            endif

        ! Store optical depth
        outod(i,j,l,ibnd) = outod(i,j,l,ibnd) + optdep

      enddo  ! non-dust species


      ! Dust optical depth
      do ispec = 1, ndust

        ! Dust index
        idst = nspecs - 1 + ispec

        ! Get Q and Reff
        ! Convert Reff um --> cm
        qext = q_dext(ispec,ibnd)
        reff = r_dust(ispec,ibnd) * 1e-4

        ! Compute box optical depth
        optdep = 0.75 * qext * dz*conc_gcm3(idst+2) / &
                 ( reff * rho(idst) )

        ! Store optical depth
        outod(i,j,l,ibnd) = outod(i,j,l,ibnd) + optdep
        ! Fine mode
!        if (ispec<=4) outod(i,j,l,nspecs+1) = outod(i,j,l,nspecs+1) + optdep

        ! Number concentration [#/cm3]
        numconc = conc_gcm3(idst+2) / &
          (rho(idst) * v_dave(ispec,ibnd)*1e-12)
        ! Extinction [km-1]
        k_ext = c_dext(ispec,ibnd) * numconc * 1e-3
        kappa_ext = kappa_ext + k_ext
        ! Scattering [km-1]
        k_sca = c_dsca(ispec,ibnd) * numconc * 1e-3
        kappa_sca = kappa_sca + k_sca
        ! Single Scattering Albedo
        totssa = totssa + k_ext * dssalb(ispec,ibnd)
        ! Asymmetry parameter
        totg = totg + k_sca * dasym(ispec,ibnd)

      enddo  ! dust species

      ! Store total single scattering albedo
      if (kappa_ext .ne. 0) then
         outssa(i,j,l,ibnd) = totssa / kappa_ext
      endif

      ! Store total asymmetry parameter
      if (kappa_sca .ne. 0) then
         outg(i,j,l,ibnd) = totg / kappa_sca
      endif

    enddo

  enddo
  enddo
  enddo

end subroutine calc_aod

!======================================================================
! ALLOC_OPTICS: Allocate OPTICS arrays
!======================================================================
subroutine alloc_optics

  allocate(rhbins(nrh))
  allocate(wlbins(nwl))
  allocate(specmr(nwl,nrh,nspecs))
  allocate(specmi(nwl,nrh,nspecs))
  allocate(mr(nrh,nspecs,nbnd))
  allocate(mi(nrh,nspecs,nbnd))
  allocate(q_ext(nrh,nspecs-1,nbnd))
  allocate(c_ext(nrh,nspecs-1,nbnd))
  allocate(c_sca(nrh,nspecs-1,nbnd))
  allocate(v_ave(nrh,nspecs-1,nbnd))
  allocate(ssalb(nrh,nspecs-1,nbnd))
  allocate(asym(nrh,nspecs-1,nbnd))
  allocate(r_eff(nrh,nspecs-1,nbnd))
  allocate(q_dext(ndust,nbnd))
  allocate(c_dext(ndust,nbnd))
  allocate(c_dsca(ndust,nbnd))
  allocate(v_dave(ndust,nbnd))
  allocate(dssalb(ndust,nbnd))
  allocate(dasym(ndust,nbnd))
  allocate(r_dust(ndust,nbnd))

end subroutine alloc_optics

!======================================================================
! ALLOC_DIST: Allocate size distributions arrays
!======================================================================
subroutine alloc_dist

  allocate(radius(nrh,nsizspc))
  allocate(rmin(nrh,nsizspc))
  allocate(rmax(nrh,nsizspc))

end subroutine alloc_dist

!======================================================================
! ALLOC_OTHER: Allocate other arrays
!======================================================================
subroutine alloc_other

  allocate(outod(ni,nj,nl,nbnd))
  allocate(outssa(ni,nj,nl,nbnd))
  allocate(outg(ni,nj,nl,nbnd))

end subroutine alloc_other

!======================================================================

end module flexaod

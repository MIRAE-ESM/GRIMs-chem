#include <define.h>
!=====================================================================
!
!     MODULE FLEXAOD
!
!        This module calculates the optical properties of aerosols
!        such as AOD, SSA, ASM based upon MIE algorithm.
!
!                                               written by G. Curci
!                                           modified by Seungun Lee
!
!=====================================================================

module flexaod

  use dao_mod,       only : airvol, bxheight, rh
  use directory_mod, only : data_dir
  use tracer_mod,    only : stt
  use tracerid_mod,  only : idtso4, idtnh4, idtnit
  use tracerid_mod,  only : idtbcpi, idtocpi, idtbcpo, idtocpo
  use tracerid_mod,  only : idtsala, idtsalc
  use tracerid_mod,  only : idtdst1, idtdst2, idtdst3, idtdst4
  use cmn_size_mod,  only : ni => iipar, nj => jjpar, nl => llpar
  use cmn_size_mod,  only : master, lonloop, latloop

  implicit none

  ! I/O units
  integer,parameter :: iou=10          ! Generic use

  ! Parameters
  integer,parameter :: nwl = 61        ! number input wavelengths
  integer,parameter :: nrh = 5         ! number of input RH bins
  integer,parameter :: nsinyuk = 109
  integer,parameter :: nbnd = 8
  integer,parameter :: nspecs = 6      ! Number of GEOS-Chem species
  integer,parameter :: ndust = 7       ! Number of GEOS-Chem DUST species
  integer,parameter :: npcoef = 100    ! Number of Legendre expansion coeff.
  integer,parameter :: nsizspc = nspecs + ndust - 1
  character(len=4),dimension(nspecs),parameter :: &
    species = (/ 'sulf', 'oc', 'bc', 'ssa', 'ssc', 'dust' /)
  real,dimension(6), parameter :: &
    hgf_bc_chin = (/ 1.0, 1.0, 1.0, 1.2, 1.4, 1.5 /)

  !--------------
  ! Optics vars
  !--------------
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
  real,dimension(:,:,:),allocatable   :: relh, boxh, aird
  real,dimension(:,:,:,:),allocatable :: outod
  real,dimension(:,:,:,:),allocatable :: outssa, outg
  real,dimension(:,:,:,:),allocatable :: k_ext, k_sca
  real,dimension(:,:,:,:),allocatable :: ssa, numconc, ga
  real,dimension(nsizspc)             :: rho, sigma, gam_a, gam_b
  real,dimension(:,:),allocatable     :: radius, rmin, rmax
  real,dimension(:,:,:,:),allocatable :: conc

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

  subroutine mie_calc

  ! Wavelenght in microns (um)
  wl_um(1) = 0.255
  wl_um(2) = 0.2925
  wl_um(3) = 0.3125
  wl_um(4) = 0.5075
  wl_um(5) = 3.135
  wl_um(6) = 1.745
  wl_um(7) = 0.96
  wl_um(8) = 2.35

  !--------------------------------------------------------------------
  ! Read Optical properties
  !--------------------------------------------------------------------

  ! Message
  if (master) then
  write(6,'("  o Reading optical properties ... ")')
  endif

  ! Allocate arrays
  call alloc_optics

  !====================================================================
  ! Switch optics source
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

          if (master) then
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

#ifdef MP
          call mpbcastr(wlbins,nwl)
          call mpbcastr(specmr,nwl*nrh*nspecs)
          call mpbcastr(specmi,nwl*nrh*nspecs)
#endif

          ! Exit immediately if BC and dust (dry data only)
          if (species(ispec)=='bc' .or. species(ispec)=='dust') exit

        enddo
      enddo

      ! Interpolation of refractive indices onto output wavelength
      do ibnd = 1, nbnd
        do ispec = 1, nspecs
          do irh = 1, nrh
            mr(irh,ispec,ibnd) = linterp( wlbins, specmr(:,irh,ispec), wl_um(ibnd) )
            mi(irh,ispec,ibnd) = linterp( wlbins, specmi(:,irh,ispec), wl_um(ibnd) )
          enddo
        enddo
      enddo

  !--------------------------------------------------------------------
  ! Overwrite dust refractive indices for STD case with
  ! those provided by Sinyuk et al. (2003)

    ! Reallocate buffer arrays
    if (allocated(wlbins)) deallocate(wlbins)
    if (allocated(specmr)) deallocate(specmr)
    if (allocated(specmi)) deallocate(specmi)
    allocate(wlbins(nsinyuk))
    allocate(specmr(nsinyuk,1,1))
    allocate(specmi(nsinyuk,1,1))

    if (master) then
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
    call mpbcastr(wlbins,nsinyuk)
    call mpbcastr(specmr,nsinyuk)
    call mpbcastr(specmi,nsinyuk)
#endif

    do ibnd = 1, nbnd
      ! Interpolation of refractive indices onto output wavelength
      mr(1,nspecs,ibnd) = linterp( wlbins, specmr(:,1,1), max(wl_um(ibnd),0.3) )
      ! Change sign of imaginary index for consistency with Mie calculations
      mi(1,nspecs,ibnd) = - linterp( wlbins, specmi(:,1,1), max(wl_um(ibnd),0.3) )
    enddo

  !--------------------------------------------------------------------

  !====================================================================

  !--------------------------------------------------------------------
  ! Read Water optical properties
  !--------------------------------------------------------------------

      mwatr = 1.33
      mwati = 0.0

  !--------------------------------------------------------------------
  ! Read size distributions
  !--------------------------------------------------------------------

  ! Message
  if (master) then
  write(6,'("  o Reading size distributions ... ")')
  endif

  ! Allocate arrays
  call alloc_dist

  if (master) then
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
  call mpbcastr(rho,nsizspc)
  call mpbcastr(radius,nrh*nsizspc)
  call mpbcastr(sigma,nsizspc)
  call mpbcastr(gam_a,nsizspc)
  call mpbcastr(gam_b,nsizspc)
  call mpbcastr(rmin,nrh*nsizspc)
  call mpbcastr(rmax,nrh*nsizspc)
#endif

  !====================================================================
  ! Hygroscopic growth

      if (master) then
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
  call mpbcastr(radius,nrh*nsizspc)
  call mpbcastr(rmin,nrh*nsizspc)
  call mpbcastr(rmax,nrh*nsizspc)
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

    if (master) then
    filename = trim(data_dir)//'/flexaod/'//'mie_tables'
    open(111, file=filename, form='unformatted', status='old', iostat=ios)
    endif

#ifdef MP
    call mpbcasti(ios,1)
#endif

    if (ios .ne. 0) then
      if (master) then
      write(6,'("  o Reading Mie tables failed !!! ")')
      endif
      go to 3000
    endif

    if (master) then
    read(111) q_ext, c_ext, c_sca, r_eff, v_ave, ssalb, asym, &
              q_dext, c_dext, c_dsca, r_dust, v_dave, dssalb, dasym
    close(111)
    endif

#ifdef MP
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
#endif

    if (master) then
    write(6,'("  o Reading Mie tables done !!! ")')
    endif

    return

3000 continue

    ! Message
    if (master) then
    write(6,'("  o Calculating Mie tables ... ")')
    endif

    ! Mie calculations
    call mie_tab

    ! Message
    if (master) then
    write(6,'("  o Mie tables done !!! ")')
    endif

    if (master) then
    filename = trim(data_dir)//'/flexaod/'//'mie_tables'
    open(111, file=filename, form='unformatted', status='unknown')
    write(111) q_ext, c_ext, c_sca, r_eff, v_ave, ssalb, asym, &
               q_dext, c_dext, c_dsca, r_dust, v_dave, dssalb, dasym
    close(111)
    endif

    ! Message
    if (master) then
    write(6,'("  o Writing Mie tables done !!! ")')
    endif

  end subroutine mie_calc

!======================================================================

  subroutine calc_aod

  if (.not. allocated(conc)) then
    call alloc_other
    outod = 0d0
    outssa = 0d0
    outg = 0d0
    conc = 0d0
  endif

  ! Tracer concentrations (kg)
  if (idtso4 .ne. 0 .and. idtnh4 .ne. 0 .and. idtnit .ne. 0) then
    conc(:,:,:,1) = stt(:,:,:,idtso4) + stt(:,:,:,idtnh4) + stt(:,:,:,idtnit)
  endif
  if (idtocpi .ne. 0) then
    conc(:,:,:,2) = stt(:,:,:,idtocpi)
  endif
  if (idtbcpi .ne. 0) then
    conc(:,:,:,3) = stt(:,:,:,idtbcpi)
  endif
  if (idtsala .ne. 0) then
    conc(:,:,:,4) = stt(:,:,:,idtsala)
  endif
  if (idtsalc .ne. 0) then
    conc(:,:,:,5) = stt(:,:,:,idtsalc)
  endif
  if (idtocpo .ne. 0) then
    conc(:,:,:,6) = stt(:,:,:,idtocpo)
  endif
  if (idtbcpo .ne. 0) then
    conc(:,:,:,7) = stt(:,:,:,idtbcpo)
  endif
  if (idtdst1 .ne. 0) then
    conc(:,:,:,8) = stt(:,:,:,idtdst1)*0.25
    conc(:,:,:,9) = stt(:,:,:,idtdst1)*0.25
    conc(:,:,:,10) = stt(:,:,:,idtdst1)*0.25
    conc(:,:,:,11) = stt(:,:,:,idtdst1)*0.25
  endif
  if (idtdst2 .ne. 0) then
    conc(:,:,:,12) = stt(:,:,:,idtdst2)
  endif
  if (idtdst3 .ne. 0) then
    conc(:,:,:,13) = stt(:,:,:,idtdst3)
  endif
  if (idtdst4 .ne. 0) then
    conc(:,:,:,14) = stt(:,:,:,idtdst4)
  endif

  ! RH (%)
  relh = rh
  ! Air volume (m3)
  aird = airvol
  ! Box height (m)
  boxh = bxheight

  if (sum(airvol) .eq. 0) return

  call interp_aod

  end subroutine calc_aod

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
    if (master) then
    write(6,'("!!! Output point must be in the range ",F6.3," - ", &
                   F6.3," !!!")') &
           minval(bins), maxval(bins)
    endif
#ifdef MP
    call mpabort
#else
    stop
#endif
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
  if (master) then
  write(6,'("   > For Wavelenth = ",F6.4," um ")') wl_um(ibnd)
  endif

  ! Non-dust species
  do ispec = 1, nspecs-1

    ! Message
    if (master) then
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
      do ip = 1, npcoef
!        pcoef(ip,irh,ispec,ibnd) = real(al1(ip)) / (2.*(real(ip)-1.)+1.)
      enddo

      ! debug
      !write(*,'(I1,2X,F6.4,2X,F5.3,2x,F6.4,1X,F5.3)')irh,q_ext(irh,ispec),r_eff(irh,ispec),ssalb(irh,ispec),al1(2)
    enddo
  enddo

  ! Dust species
  ! Message
  if (master) then
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
    do ip = 1, npcoef
!      dpcoef(ip,ispec,ibnd) = real(al1(ip)) / (2.*(real(ip)-1.)+1.)
    enddo

    ! debug
    !write(*,'(I1,2X,F6.4,2X,F5.3,2x,F5.3,8(1X,F5.3))')irh,q_dext(ispec),r_dust(ispec),dssalb(ispec),al1(1:8)
  enddo

  enddo

end subroutine mie_tab

!======================================================================
! INTERP_AOD: Calculate RH-dependent AOD from Mie table
!======================================================================
subroutine interp_aod

  integer :: i, j, l, ispec, irh, iphob, idst, r, ip
  integer :: n, icount
  real :: rh, dz, fwet, weight
  real :: reff, qext, optdep, g, ss
  real :: scaleq, scaler, scaleod
  real :: x, y
  real :: kappa_ext, kappa_sca, totssa, totg
  real,dimension(npcoef) :: pci
  real,dimension(nsizspc+2) :: conc_gcm3
  real,dimension(nrh) :: rw, qw, gw, ssw
  real,dimension(npcoef,nrh) :: pcw

  ! Loop over horizontal grid
  outod = 0.0
  outssa = 0.0
  outg = 0.0
  do ibnd = 1,nbnd
  do j = 1,latloop
  do i = 1,lonloop(j)

    ! Loop over levels
    do l = 1,nl

      ! Aerosol concentration (g/cm3)
      conc_gcm3 = conc(i,j,l,:) / aird(i,j,l) * 1e-3

      ! Relative Humidity (%)
      rh = relh(i,j,l)

      ! Box Height (m --> cm)
      dz = boxh(i,j,l) * 1e2


      ! Reset total variables
      kappa_ext = 0.0
      kappa_sca = 0.0
      totssa = 0.0
      totg = 0.0


        ! Select RH bin
        do irh = 2,nrh
          if (rh<rhbins(irh)) exit
        enddo
        irh = irh - 1

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

              ! Wet phase function coeffs
              do ip = 1, npcoef
                !pcw(ip,r) = pcoef(ip,r,ispec)*fwet + &
                !            pcoef(ip,1,ispec)*(1.d0-fwet)
!                pcw(ip,r) = pcoef(ip,r,ispec,ibnd)
              enddo

            enddo  ! RH bins

            ! Interpolate optical parameters
            if (irh<nrh) then

              ! Interpolation weight
              weight = (rh-rhbins(irh)) / (rhbins(irh+1)-rhbins(irh))
              if ( weight > 1.0d0 ) weight = 1.0d0

              ! Interpolate radius and Q scaling factor
              scaleq = (weight*qw(irh+1) + (1.d0-weight)*qw(irh)) / qw(1)
              reff   =  weight*rw(irh+1) + (1.d0-weight)*rw(irh)

              ! Interpolate single-scattering albedo
              ss =  weight*ssw(irh+1) + (1.d0-weight)*ssw(irh)

              ! Interpolate asymmetry parameter
              g =  weight*gw(irh+1) + (1.d0-weight)*gw(irh)

              ! Interpolate phase function coeffs
              do ip = 1, npcoef
!                pci(ip) = weight*pcw(ip,irh+1) + (1.d0-weight)*pcw(ip,irh)
              enddo

            ! Last RH bin: do not interpolate
            else

              ! Radius and Q scaling factor
              scaleq = qw(nrh) / qw(1)
              reff = rw(nrh)

              ! Single-scattering albedo
              ss =  ssw(irh)

              ! Asymmetry parameter
              g =  gw(irh)

              ! Phase function coeffs
              do ip = 1, npcoef
!                pci(ip) = pcw(ip,irh)
              enddo

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
            numconc(ispec,i,j,l) = conc_gcm3(ispec) / &
              (rho(ispec) * v_ave(1,ispec,ibnd)*1e-12)
            ! Extinction coefficient [km-1]
            !  Cext = extinction cross-section (um^2)
            !  Ncon = number concentration (cm^-3)
            !  Kext = extinction coefficient (km^-1)
            k_ext(ispec,i,j,l) = c_ext(1,ispec,ibnd) * &
              numconc(ispec,i,j,l) * scaleod*1e-3
            kappa_ext = kappa_ext + k_ext(ispec,i,j,l)
            ! Scattering coefficient [km-1]
            k_sca(ispec,i,j,l) = c_sca(1,ispec,ibnd) * &
              numconc(ispec,i,j,l) * scaleod*1e-3
            kappa_sca = kappa_sca + k_sca(ispec,i,j,l)
            ! Back-Scattering coefficient [km-1 srad-1]
!            k_bac(ispec,i,j,l) = c_bac(1,ispec,ibnd) * &
!              numconc(ispec,i,j,l) * scaleod*1e-3
            ! Single Scattering Albedo
            ssa(ispec,i,j,l) = ss
            totssa = totssa + k_ext(ispec,i,j,l) * ss
            ! Asymmetry factor
            ga(ispec,i,j,l) = g
            totg = totg + k_sca(ispec,i,j,l) * g
            ! Phase function coefficients
            do ip = 1, npcoef
!              pc(ip,ispec,i,j,l) = pci(ip)
            enddo

            ! Hydrophobic species
            if (iphob>0) then
              ! Number concentration [#/cm3]
              numconc(iphob,i,j,l) = conc_gcm3(iphob) / &
                (rho(ispec) * v_ave(1,ispec,ibnd)*1e-12)
              ! Extinction [km-1]
              k_ext(iphob,i,j,l) = c_ext(1,ispec,ibnd) * &
                numconc(iphob,i,j,l) * 1e-3
              kappa_ext = kappa_ext + k_ext(iphob,i,j,l)
              ! Scattering [km-1]
              k_sca(iphob,i,j,l) = c_sca(1,ispec,ibnd) * &
                numconc(iphob,i,j,l) * 1e-3
              kappa_sca = kappa_sca + k_sca(iphob,i,j,l)
              ! Back-Scattering [km-1 srad-1]
!              k_bac(iphob,i,j,l) = c_bac(1,ispec,ibnd) * &
!                numconc(iphob,i,j,l) * 1e-3
              ! Single Scattering Albedo
              ssa(iphob,i,j,l) = ssalb(1,ispec,ibnd)
              totssa = totssa + k_ext(iphob,i,j,l) * ssalb(1,ispec,ibnd)
              ! Asymmetry factor
              ga(iphob,i,j,l) = asym(1,ispec,ibnd)
              totg = totg + k_sca(iphob,i,j,l) * asym(1,ispec,ibnd)
              ! Phase function coefficients
              do ip = 1, npcoef
!                pc(ip,iphob,i,j,l) = pcoef(ip,1,ispec,ibnd)
              enddo
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
        numconc(idst+2,i,j,l) = conc_gcm3(idst+2) / &
          (rho(idst) * v_dave(ispec,ibnd)*1e-12)
        ! Extinction [km-1]
        k_ext(idst+2,i,j,l) = c_dext(ispec,ibnd) * numconc(idst+2,i,j,l) * 1e-3
        kappa_ext = kappa_ext + k_ext(idst+2,i,j,l)
        ! Scattering [km-1]
        k_sca(idst+2,i,j,l) = c_dsca(ispec,ibnd) * numconc(idst+2,i,j,l) * 1e-3
        kappa_sca = kappa_sca + k_sca(idst+2,i,j,l)
        ! Back-Scattering [km-1 srad-1]
!        k_bac(idst+2,i,j,l) = c_dbac(ispec,ibnd) * numconc(idst+2,i,j,l) * 1e-3
        ! Single Scattering Albedo
        ssa(idst+2,i,j,l) = dssalb(ispec,ibnd)
        totssa = totssa + k_ext(idst+2,i,j,l) * dssalb(ispec,ibnd)
        ! Asymmetry parameter
        ga(idst+2,i,j,l) = dasym(ispec,ibnd)
        totg = totg + k_sca(idst+2,i,j,l) * dasym(ispec,ibnd)
        ! Phase function coefficients
        do ip = 1, npcoef
!          pc(ip,idst+2,i,j,l) = dpcoef(ip,ispec,ibnd)
        enddo

      enddo  ! dust species

      ! Store total single scattering albedo
      outssa(i,j,l,ibnd) = totssa / kappa_ext

      ! Store total asymmetry parameter
      outg(i,j,l,ibnd) = totg / kappa_sca

    enddo  ! levels

  enddo
  enddo
  enddo

end subroutine interp_aod

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

  allocate(conc(ni,nj,nl,nsizspc+2))
  allocate(relh(ni,nj,nl))
  allocate(aird(ni,nj,nl))
  allocate(boxh(ni,nj,nl))
  allocate(outod(ni,nj,nl,nbnd))
  allocate(outssa(ni,nj,nl,nbnd))
  allocate(outg(ni,nj,nl,nbnd))
  allocate(k_ext(nsizspc+2,ni,nj,nl))
  allocate(k_sca(nsizspc+2,ni,nj,nl))
  allocate(ssa(nsizspc+2,ni,nj,nl))
  allocate(ga(nsizspc+2,ni,nj,nl))
  allocate(numconc(nsizspc+2,ni,nj,nl))

end subroutine alloc_other

end module flexaod

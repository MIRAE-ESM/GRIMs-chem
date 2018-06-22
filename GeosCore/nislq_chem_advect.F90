#include <define.h>
   subroutine nislq_chem_advect(deltim,pt,ut,vt,pdot)
#ifndef RMP
!-------------------------------------------------------------------------------
!
! a routine to do non-iteration semi-Lagrangain finite volume advection
! considering mass advection together with gases
! contact: hann-ming henry juang
!
! <in : top to bottom>
! deltim  time step from n to n+1 divided by 2
! pt      surface pressure in cb
! ut      horizontal u wind in m/s
! vt      horizontal v wind in m/s
! pdot    vertical wind in dp/dt in cb
!
!-------------------------------------------------------------------------------
   use constant   , only : rrerth_
   use comfgrid   , only : rbs2
   use comfver    , only : ak5,bk5
   use paramodel  , only : LONF2S,LATG2S,lonf2_,latg2_,levs_,levh_,LEVSS
   use comfspec_vr, only : qm,z
   use comio      , only : iope
#ifdef MP
   use commpi     , only : nsize=>nrow,ncol,mype,latlen,latdef,latstr
   use paramodel  , only : lonf2p_,latg2p_,levsp_
#else
   use comfgrid   , only : latdef
#endif
   use nislq      , only : nx,lev,my,my_max,ncld,nlevs,nlevsp                 ,&
                           lonfull,latfull,lonpart,latpart,mylonlen           ,&
                           cyclic_cell_intpx                                  ,&
                           cyclic_cell_massadvx                               ,&
                           cyclic_cell_massadvy                               ,&
                           vertical_cell_advect
   use tracer_mod , only : nt=>n_tracers, stt
   use cmn_size_mod
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   real, intent(in)                                        ::  deltim
   real, intent(in)   , dimension(LONF2S,        LATG2S)   ::  pt
   real, intent(in)   , dimension(LONF2S,levs_  ,LATG2S)   ::  ut,vt
   real, intent(in)   , dimension(LONF2S,levs_+1,LATG2S)   ::  pdot
!
! local variables
!
   integer            , parameter                          ::  mass=0
#ifndef MP
   integer            , parameter                          ::  nsize=1
#endif
   integer                                                 ::  jjend,j,j1,j2  ,&
                                                               jj,lat         ,&
                                                               i,lonsd        ,&
                                                               k,kk,t
   real               , dimension(LONF2S ,levs_+1       )     ::  ppi,pdot2
   real               , dimension(LONF2S ,nt*levs_ ,LATG2S )  ::  qp
   real               , dimension(lonf2_ ,nt*LEVSS ,LATG2S )  ::  qt
   real               , dimension(LONF2S ,levs_ ,LATG2S )     ::  qp1
   real               , dimension(lonf2_ ,LEVSS ,LATG2S )     ::  qt1
#ifdef MP
   real               , dimension(lonf2_ ,levsp_,latg2p_)     ::  utp,vtp
   real               , dimension(lonf2_ ,nt*levsp_,latg2p_)  ::  qpp
   real               , dimension(LONF2S ,nt*levs_ ,LATG2S )  ::  qtp
   real               , dimension(lonf2_ ,levsp_,latg2p_)     ::  qpp1
   real               , dimension(LONF2S ,levs_ ,LATG2S )     ::  qtp1
#endif
   real               , dimension(LONF2S ,levs_            )  ::  qtn
   !
   real               , dimension(lonfull,LEVSS ,latpart)     ::  uulon,vvlon
   real               , dimension(latfull,LEVSS ,lonpart)     ::  vvlat
   real               , dimension(lonfull,nt*LEVSS ,latpart)  ::  qqlon,rrlon
   real               , dimension(latfull,nt*LEVSS ,lonpart)  ::  qqlat,rrlat
   real               , dimension(lonfull,LEVSS ,latpart)     ::  qqlon1,rrlon1
   real               , dimension(latfull,LEVSS ,lonpart)     ::  qqlat1,rrlat1
!
! initialize
!
   qp=0. ;  qt=0. ;  ppi=0.
#ifdef MP
   utp=0.;  vtp=0.;  qpp=0.;  qtp=0.
#endif
   uulon=0. ;  vvlon=0. ;  vvlat=0.
   qqlon=0. ;  rrlon=0. ;  qqlat=0.  ; rrlat=0.
!
! latitude band
!
#ifdef MP
   jjend=latlen(mype)
#else
   jjend=latg2_
#endif
!
! qp[t2b], stt=[b2t]
!
!$omp parallel
!$omp do collapse(2) private(t,k,j,i)
   do t = 1,nt
     do k = 1,llpar
       do j = 1,jjpar
         do i = 1,iipar
           qp(i,llpar*t+1-k,j) = stt(i,j,k,t)
         enddo
       enddo
     enddo
   enddo
#ifdef MP
!
! transpose z-full to x-full
!
!$omp master
   call mpnx2nk(ut,lonf2p_,levs_,utp,lonf2_,levsp_,latg2p_,levs_,levsp_,       &
                                                                 1,1,1)
   call mpnx2nk(vt,lonf2p_,levs_,vtp,lonf2_,levsp_,latg2p_,levs_,levsp_,       &
                                                                 1,1,1)
   do t = 1,nt
     kk = levs_*(t-1)
     do j = 1,latg2p_
       do k = 1,levs_
         do i = 1,lonf2p_
           qp1(i,k,j) = qp(i,kk+k,j)
         enddo
       enddo
     enddo
     call mpnx2nk(qp1,lonf2p_,levs_,qpp1,lonf2_,levsp_,latg2p_,levs_,levsp_,   &
                                                                 1,1,1)
     kk = LEVSS*(t-1)
     do j = 1,latg2p_
       do k = 1,LEVSS
         do i = 1,lonf2_
           qpp(i,kk+k,j) = qpp1(i,k,j)
         enddo
       enddo
     enddo
   enddo
!$omp end master
!$omp barrier
#define UU utp
#define VV vtp
#define QQ qpp
#else
#define UU ut
#define VV vt
#define QQ qp
#endif /* MP end */
!
! first mass conserving interpolation from reduced grid to full grid
!
!$omp do private(jj,j1,j2,lat,lonsd,k,i)
   do jj = 1,jjend
     j1=2*jj-1     ! N.H.
     j2=2*jj       ! S.H.
     lat=latdef(jj)
     lonsd=nx
     ! u,v at time step n (convert dynamics to SL grid)
     do k = 1,LEVSS
       do i = 1,lonsd
         uulon(i,k,j1) = UU(i      ,k,jj) * (rrerth_*sqrt(rbs2(jj)))
         uulon(i,k,j2) = UU(lonsd+i,k,jj) * (rrerth_*sqrt(rbs2(jj)))
         vvlon(i,k,j1) = VV(i      ,k,jj) * (rrerth_)
         vvlon(i,k,j2) = VV(lonsd+i,k,jj) * (rrerth_)
       enddo
     enddo
   enddo
! ---------------------------------------------------------------------
! mpi para from horizontal full grid to meridional full grid
! ---------------------------------------------------------------------
!
! para vvlon to vvlat
!
!$omp master
   call nislq_transpose_we2ns(vvlon,vvlat,LEVSS,nsize)
!$omp end master
!$omp barrier
!
! first mass conserving interpolation from reduced grid to full grid
!
!$omp do collapse(2) private(t,jj,kk,j1,j2,lat,lonsd,k,i)
 do t = 1,nt
   do jj = 1,jjend
     kk = LEVSS*(t-1)
     j1=2*jj-1     ! N.H.
     j2=2*jj       ! S.H.
     lat=latdef(jj)
     lonsd=nx
     ! at time step n-1 (convert dynamics to SL grid)
     do k = 1,LEVSS
       do i = 1,lonsd
         qqlon(i,kk+k,j1) = QQ(i      ,kk+k,jj)
         qqlon(i,kk+k,j2) = QQ(lonsd+i,kk+k,jj)
       enddo
     enddo
#undef UU
#undef VV
#undef QQ
!
     do k = 1,LEVSS
       do i = 1,lonfull
         rrlon(i,kk+k,j1) = qqlon(i,kk+k,j1)
         rrlon(i,kk+k,j2) = qqlon(i,kk+k,j2)
       enddo
     enddo
!
! first set positive advection in horziontal direction with mass conserving
!
     call cyclic_cell_massadvx(LEVSS,1,deltim,                                 &
                                          uulon(1,1,j1),rrlon(1,kk+1,j1),mass)
     call cyclic_cell_massadvx(LEVSS,1,deltim,                                 &
                                          uulon(1,1,j2),rrlon(1,kk+1,j2),mass)
   enddo
 enddo
! ---------------------------------------------------------------------
! mpi para from horizontal full grid to meridional full grid
! ---------------------------------------------------------------------
!
! para qqlon and rrlon to qqlat and rrlat
!
!$omp master
   do t = 1,nt
     kk = LEVSS*(t-1)
     do j = 1,latpart
       do k = 1,LEVSS
         do i = 1,lonfull
           qqlon1(i,k,j) = qqlon(i,kk+k,j)
           rrlon1(i,k,j) = rrlon(i,kk+k,j)
         enddo
       enddo
     enddo
     call nislq_transpose_we2ns(qqlon1,qqlat1,LEVSS,nsize)
     call nislq_transpose_we2ns(rrlon1,rrlat1,LEVSS,nsize)
     do i = 1,lonpart
       do k = 1,LEVSS
         do j = 1,latfull
           qqlat(j,kk+k,i) = qqlat1(j,k,i)
           rrlat(j,kk+k,i) = rrlat1(j,k,i)
         enddo
       enddo
     enddo
   enddo
!$omp end master
!$omp barrier
! ---------------------------------------------------------------------
! ------------------- in meridional great circle ----------------------
! ---------------------------------------------------------------------
!$omp do collapse(2) private(t,i,kk)
 do t = 1,nt
   do i = 1,mylonlen
   kk = LEVSS*(t-1)
! 
! first set advection in meridional direction in great circle through two poles
!
     call cyclic_cell_massadvy(LEVSS,1,deltim,vvlat(1,1,i),rrlat(1,kk+1,i),mass)
!
! second set advection in meridional direction in great circle through two poles
!
     call cyclic_cell_massadvy(LEVSS,1,deltim,vvlat(1,1,i),qqlat(1,kk+1,i),mass)

   enddo
 enddo
!
! ----------------------------------------------------------------------
! mpi para from meridional direction to horizontal directory 
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon
!
!$omp master
   do t = 1,nt
     kk = LEVSS*(t-1)
     do i = 1,lonpart
       do k = 1,LEVSS
         do j = 1,latfull
           qqlat1(j,k,i) = qqlat(j,kk+k,i)
           rrlat1(j,k,i) = rrlat(j,kk+k,i)
         enddo
       enddo
     enddo
     call nislq_transpose_ns2we(qqlat1,qqlon1,LEVSS,nsize)
     call nislq_transpose_ns2we(rrlat1,rrlon1,LEVSS,nsize)
     do j = 1,latpart
       do k = 1,LEVSS
         do i = 1,lonfull
           qqlon(i,kk+k,j) = qqlon1(i,k,j)
           rrlon(i,kk+k,j) = rrlon1(i,k,j)
         enddo
       enddo
     enddo
   enddo
!$omp end master
!$omp barrier
! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!     print *,' nislq adv loop in x for last '
!
!$omp do collapse(2) private(t,jj,kk,j1,j2,lat,lonsd)
 do t = 1,nt
   do jj = 1,jjend
     kk = LEVSS*(t-1)
     j1=2*jj-1     ! N.H.
     j2=2*jj       ! S.H.
     lat=latdef(jj)
     lonsd=nx
!
! second set advection in x for the second of the pair
!
     call cyclic_cell_massadvx(LEVSS,1,deltim,                                 &
                                          uulon(1,1,j1),qqlon(1,kk+1,j1),mass)
     call cyclic_cell_massadvx(LEVSS,1,deltim,                                 &
                                          uulon(1,1,j2),qqlon(1,kk+1,j2),mass)
     do k = 1,LEVSS
       do i = 1,lonsd
         rrlon(i,kk+k,j1) = 0.5 * ( qqlon(i,kk+k,j1) + rrlon(i,kk+k,j1) )
         rrlon(i,kk+k,j2) = 0.5 * ( qqlon(i,kk+k,j2) + rrlon(i,kk+k,j2) )
       enddo
     enddo
!
! convert SL to dynamics grid
!
     do k = 1,LEVSS
       do i = 1,lonsd
         qt(i      ,kk+k,jj)=rrlon(i,kk+k,j1)
         qt(lonsd+i,kk+k,jj)=rrlon(i,kk+k,j2)
       enddo
     enddo
   enddo
 enddo
#ifdef MP
!
! transpose x-full to z-full
!
!$omp master
   do t = 1,nt
     kk = LEVSS*(t-1)
     do j = 1,latg2p_
       do k = 1,LEVSS
         do i = 1,lonf2_
           qt1(i,k,j) = qt(i,kk+k,j)
         enddo
       enddo
     enddo
     call mpnk2nx(qt1,lonf2_,levsp_,qtp1,lonf2p_,levs_,latg2p_,levsp_,levs_,   &
                                                                  1,1,1)
     kk = levs_*(t-1)
     do j = 1,latg2p_
       do k = 1,levs_
         do i = 1,lonf2p_
           qtp(i,kk+k,j) = qtp1(i,k,j)
         enddo
       enddo
     enddo
   enddo
!$omp end master
!$omp barrier
#define QT qtp
#else
#define QT qt
#endif /* MP end */
! --------------------------------------------------------------
! ----------- compute vertical advection and total ------------
! --------------------------------------------------------------
!$omp do collapse(2) private(t,j,kk,lonsd,k,i,ppi,pdot2,qtn)
 do t = 1,nt
   do j = 1,jjend
     kk = levs_*(t-1)
     lonsd=LONF2S
!
!    pressure (top to bottom)
!
     do k = 1,levs_+1
       do i = 1,lonsd
         ppi(i,k)=ak5(k)+bk5(k)*pt(i,j)
         pdot2(i,k)=pdot(i,k,j)
       enddo
     enddo
!
     do k = 1,levs_
       do i = 1,lonsd
         qtn(i,k)=QT(i,kk+k,j)
       enddo
     enddo
#undef QT
!
! vertical advection with mass conserving positive advection
!
     call vertical_cell_advect(lonsd,LONF2S,lev,1,deltim,ppi,pdot2,qtn,mass)
!
!    q update at time step n+1 (bottom to top)
!
     do k = 1,levs_
       do i = 1,iipar
         stt(i,j,levs_+1-k,t)=qtn(i,k)
       enddo
     enddo
   enddo
 enddo
!$omp end parallel
!
!
#endif
   return
   end subroutine nislq_chem_advect

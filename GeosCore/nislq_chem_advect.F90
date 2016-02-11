#include <define.h>
   subroutine nislq_chem_advect(deltim,pt,ut,vt,pdot,q1,q3)
!-------------------------------------------------------------------------------
!
! a routine to do non-iteration semi-Lagrangain finite volume advection
! considering mass advection together with gases
! contact: hann-ming henry juang
!
! <in : top to bottom>
! deltim  time step from n to n+1
! pt      surface pressure in cb
! ut      horizontal u wind scaled without earth radius cos(phi)^2d(lamda)/dt
! vt      horizontal v wind scaled without earth radius cos(phi)d(phi)/dt
! pdot    vertical wind in dp/dt in cb
!
! <in : bottom to top>
! q1      n-1 time step q
!
! <out: bottom to top>
! q3      n+1 time step q (bottom to top)
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
   use tracer_mod , only : ntr=>n_tracers
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
   real, intent(in)   , dimension(LONF2S,levs_*ntr  ,LATG2S)   ::  q1
   real, intent(out)  , dimension(LONF2S,levs_*ntr  ,LATG2S)   ::  q3
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
                                                               k,kk
   real               , dimension(LONF2S ,levs_+1       )  ::  ppi,pdot2
   real               , dimension(LONF2S ,levs_*ntr ,LATG2S )  ::  qp
   real               , dimension(lonf2_ ,LEVSS*ntr,LATG2S )  ::  qt
#ifdef MP
   real               , dimension(lonf2_ ,levsp_,latg2p_)  ::  utp,vtp
   real               , dimension(lonf2_ ,LEVSS*ntr,latg2p_)  ::  qpp
   real               , dimension(LONF2S ,levs_*ntr ,LATG2S )  ::  qtp
#endif
   real               , dimension(LONF2S ,levs_*ntr         )  ::  qtn
   !
   real               , dimension(lonfull,LEVSS ,latpart)  ::  uulon,vvlon
   real               , dimension(latfull,LEVSS ,lonpart)  ::  vvlat
   real               , dimension(lonfull,LEVSS*ntr,latpart)  ::  qqlon,rrlon
   real               , dimension(latfull,LEVSS*ntr,lonpart)  ::  qqlat,rrlat
!
! initialize
!
   q3=0.
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
! qp[t2b], q1=[b2t]
!
   forall(i=1:LONF2S,k=1:levs_*ntr,j=1:LATG2S) qp(i,k,j)=q1(i,levs_*ntr+1-k,j)
#ifdef MP
!
! transpose z-full to x-full
!
   call mpnx2nk(ut,lonf2p_,levs_,utp,lonf2_,levsp_,latg2p_,levs_,levsp_,       &
                                                                 1,1,1)
   call mpnx2nk(vt,lonf2p_,levs_,vtp,lonf2_,levsp_,latg2p_,levs_,levsp_,       &
                                                                 1,1,1)
   call mpnx2nk(qp,lonf2p_,levs_*ntr,qpp,lonf2_,LEVSS*ntr,latg2p_,levs_,levsp_,       &
                                                                 1,1,ntr)
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
   do jj = 1,jjend
     j1=2*jj-1     ! N.H.
     j2=2*jj       ! S.H.
     lat=latdef(jj)
     lonsd=nx
     ! u,v at time step n (convert dynamics to SL grid)
     do k = 1,LEVSS
       do i = 1,lonsd
         uulon(i,k,j1) = UU(i      ,k,jj) * (rrerth_*rbs2(jj))
         uulon(i,k,j2) = UU(lonsd+i,k,jj) * (rrerth_*rbs2(jj))
         vvlon(i,k,j1) = VV(i      ,k,jj) * (rrerth_*sqrt(rbs2(jj)))
         vvlon(i,k,j2) = VV(lonsd+i,k,jj) * (rrerth_*sqrt(rbs2(jj)))
       enddo
     enddo
     ! at time step n-1 (convert dynamics to SL grid)
     do k = 1,LEVSS*ntr
       do i = 1,lonsd
         qqlon(i,k,j1) = QQ(i      ,k,jj)
         qqlon(i,k,j2) = QQ(lonsd+i,k,jj)
       enddo
     enddo
#undef UU
#undef VV
#undef QQ
!
     do k = 1,LEVSS*ntr
       do i = 1,lonfull
         rrlon(i,k,j1) = qqlon(i,k,j1)
         rrlon(i,k,j2) = qqlon(i,k,j2)
       enddo
     enddo
!
! first set positive advection in horziontal direction with mass conserving
!
     call cyclic_cell_massadvx(LEVSS,ntr,deltim,                              &
                                          uulon(1,1,j1),rrlon(1,1,j1),mass)
     call cyclic_cell_massadvx(LEVSS,ntr,deltim,                              &
                                          uulon(1,1,j2),rrlon(1,1,j2),mass)
   enddo
! ---------------------------------------------------------------------
! mpi para from horizontal full grid to meridional full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!
   call nislq_transpose_we2ns(vvlon,vvlat,LEVSS ,nsize)
   call nislq_transpose_we2ns(qqlon,qqlat,LEVSS*ntr,nsize)
   call nislq_transpose_we2ns(rrlon,rrlat,LEVSS*ntr,nsize)
! ---------------------------------------------------------------------
! ------------------- in meridional great circle ----------------------
! ---------------------------------------------------------------------
   do i = 1,mylonlen
! 
! first set advection in meridional direction in great circle through two poles
!
     call cyclic_cell_massadvy(LEVSS,ntr,deltim,vvlat(1,1,i),rrlat(1,1,i),mass)
!
! second set advection in meridional direction in great circle through two poles
!
     call cyclic_cell_massadvy(LEVSS,ntr,deltim,vvlat(1,1,i),qqlat(1,1,i),mass)

   enddo
!
! ----------------------------------------------------------------------
! mpi para from meridional direction to horizontal directory 
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon
!
   call nislq_transpose_ns2we(qqlat,qqlon,LEVSS*ntr,nsize)
   call nislq_transpose_ns2we(rrlat,rrlon,LEVSS*ntr,nsize)
! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!     print *,' nislq adv loop in x for last '
!
   do jj = 1,jjend
     j1=2*jj-1     ! N.H.
     j2=2*jj       ! S.H.
     lat=latdef(jj)
     lonsd=nx
!
! second set advection in x for the second of the pair
!
     call cyclic_cell_massadvx(LEVSS,ntr,deltim,                              &
                                          uulon(1,1,j1),qqlon(1,1,j1),mass)
     call cyclic_cell_massadvx(LEVSS,ntr,deltim,                              &
                                          uulon(1,1,j2),qqlon(1,1,j2),mass)
     do k = 1,LEVSS*ntr
       do i = 1,lonsd
         rrlon(i,k,j1) = 0.5 * ( qqlon(i,k,j1) + rrlon(i,k,j1) )
         rrlon(i,k,j2) = 0.5 * ( qqlon(i,k,j2) + rrlon(i,k,j2) )
       enddo
     enddo
!
! convert SL to dynamics grid
!
     do k = 1,LEVSS*ntr
       do i = 1,lonsd
         qt(i      ,k,jj)=rrlon(i,k,j1)
         qt(lonsd+i,k,jj)=rrlon(i,k,j2)
       enddo
     enddo
   enddo
#ifdef MP
!
! transpose x-full to z-full
!
   call mpnk2nx(qt,lonf2_,LEVSS*ntr,qtp,lonf2p_,levs_*ntr,latg2p_,levsp_,levs_,       &
                                                                  1,1,ntr)
#define QT qtp
#else
#define QT qt
#endif /* MP end */
! --------------------------------------------------------------
! ----------- compute vertical advection and total ------------
! --------------------------------------------------------------
   do j = 1,jjend
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
     do k = 1,levs_*ntr
       do i = 1,lonsd
         qtn(i,k)=QT(i,k,j)
       enddo
     enddo
#undef QT
!
! vertical advection with mass conserving positive advection
!
     call vertical_cell_advect(lonsd,LONF2S,lev,ntr,deltim,ppi,pdot2,qtn,mass)
!
!    q update at time step n+1 (bottom to top)
!
     do k = 1,levs_*ntr
       do i = 1,lonsd
         q3(i,levs_*ntr+1-k,j)=qtn(i,k)
       enddo
     enddo
   enddo
!
!
   return
   end subroutine nislq_chem_advect

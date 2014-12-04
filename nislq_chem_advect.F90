#include <define.h>
   subroutine nislq_chem_advect(deltim,pt,ut,vt,pdot,q1,q3)
#ifdef NISLQ_MONO
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
! q1      n-1 time step q
!
! <out: bottom to top>
! q3      n+1 time step q (bottom to top)
!
!-------------------------------------------------------------------------------
   use constant   , only : rrerth_
   use comfgrid   , only : rbs2
#ifdef HYBRID
   use comfver    , only : ak5,bk5
#else
#ifdef DFS
   use dfsvar     , only : si=>sigmafull
#else
   use comfver    , only : si
#endif
#endif /* HYBRID end */
#ifdef DFS
   use dfsvar     , only : ib,jbw,iba,levs,levh,levsp,psl1,psl3,iope,latdef
#ifdef MP
   use commpi     , only : nsize=>nrow,ncol,mype,latlen,latstr
#endif
#else /* SPH */
   use paramodel  , only : LONF2S,LATG2S,lonf2_,latg2_,levs_,levh_,LEVSS
   use comfspec_vr, only : qm,z
   use comio      , only : iope
#ifdef MP
   use commpi     , only : nsize=>nrow,ncol,mype,latlen,latdef,latstr
   use paramodel  , only : lonf2p_,latg2p_,levsp_
#else
   use comfgrid   , only : latdef
#endif
#ifdef REDUCE_GRID
   use comreduce  , only : lonfd
#ifdef MP
   use comreduce  , only : lonfdp
#endif
#endif /* REDUCE_GRID end */
#endif /* DFS end */
   use nislq      , only : nx,lev,my,my_max,ncld,nlevs,nlevsp                 ,&
                           lonfull,latfull,lonpart,latpart,mylonlen           ,&
                           cyclic_cell_intpx                                  ,&
                           cyclic_cell_massadvx                               ,&
                           cyclic_cell_massadvy                               ,&
                           vertical_cell_advect
#ifdef DCMIP
   use dcmip_grims, only : ptop
#endif
#ifdef CHEM
   use tracer_mod , only : ntr=>n_tracers
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   real, intent(in)                                        ::  deltim
#ifdef DFS
   real, intent(in)   , dimension(ib,jbw)                  ::  pt
   real, intent(in)   , dimension(ib,jbw,levs)             ::  ut,vt
   real, intent(in)   , dimension(ib,jbw,levs+1)           ::  pdot
   real, intent(in)   , dimension(ib,jbw,levs*ntr)             ::  q1
   real, intent(out)  , dimension(ib,jbw,levs*ntr)             ::  q3
#define LONF2S ib
#define LATG2S jbw
#define LEVSS levsp
#define levs_ levs
#define levh_ levh
#else
   real, intent(in)   , dimension(LONF2S,        LATG2S)   ::  pt
   real, intent(in)   , dimension(LONF2S,levs_  ,LATG2S)   ::  ut,vt
   real, intent(in)   , dimension(LONF2S,levs_+1,LATG2S)   ::  pdot
   real, intent(in)   , dimension(LONF2S,levs_*ntr  ,LATG2S)   ::  q1
   real, intent(out)  , dimension(LONF2S,levs_*ntr  ,LATG2S)   ::  q3
#endif
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
#ifdef DFS
   real               , dimension(ib      , levs+1      )  ::  ppi,pdot2
   real               , dimension(ib , jbw, ntr*levs_       )  ::  qp
   real               , dimension(iba, jbw, ntr*LEVSS      )  ::  qt
#ifdef MP
   real               , dimension(iba, jbw, levsp       )  ::  utp,vtp
   real               , dimension(iba, jbw, ntr*LEVSS      )  ::  qpp
   real               , dimension(ib , jbw, ntr*levs_       )  ::  qtp
#endif
   real               , dimension(ib      , ntr*levs_       )  ::  qtn
#else /* SPH */
   real               , dimension(LONF2S ,levs_+1       )  ::  ppi,pdot2
   real               , dimension(LONF2S ,ntr*levs_ ,LATG2S )  ::  qp
   real               , dimension(lonf2_ ,ntr*LEVSS,LATG2S )  ::  qt
#ifdef MP
   real               , dimension(lonf2_ ,levsp_,latg2p_)  ::  utp,vtp
   real               , dimension(lonf2_ ,ntr*LEVSS,latg2p_)  ::  qpp
   real               , dimension(LONF2S ,ntr*levs_ ,LATG2S )  ::  qtp
#endif
   real               , dimension(LONF2S ,ntr*levs_         )  ::  qtn
#endif /* DFS end */
   !
   real               , dimension(lonfull,LEVSS ,latpart)  ::  uulon,vvlon
   real               , dimension(latfull,LEVSS ,lonpart)  ::  vvlat
   real               , dimension(lonfull,ntr*LEVSS,latpart)  ::  qqlon,rrlon
   real               , dimension(latfull,ntr*LEVSS,lonpart)  ::  qqlat,rrlat
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
! qp=q1
!
   qp=q1
#ifdef MP
!
! transpose z-full to x-full
!
#ifdef DFS
   call mpxy2yz(ut,ib,levs ,utp,iba,levsp ,jbw,levs,levsp,1)
   call mpxy2yz(vt,ib,levs ,vtp,iba,levsp ,jbw,levs,levsp,1)
   call mpxy2yz(qp,ib,ntr*levs_,qpp,iba,ntr*LEVSS,jbw,levs,levsp,ntr)
#else
   call mpnx2nk(ut,lonf2p_,levs_,utp,lonf2_,levsp_,latg2p_,levs_,levsp_,       &
                                                                 1,1,1)
   call mpnx2nk(vt,lonf2p_,levs_,vtp,lonf2_,levsp_,latg2p_,levs_,levsp_,       &
                                                                 1,1,1)
   call mpnx2nk(qp,lonf2p_,ntr*levs_,qpp,lonf2_,ntr*LEVSS,latg2p_,levs_,levsp_,       &
                                                                 1,1,ntr)
#endif
#define UU utp
#define VV vtp
#define QQ qpp
#else
#define UU ut
#define VV vt
#define QQ qp
#endif /* MP end */
#ifdef SLDBG
   if( iope ) then
     print *,' enter nislq_advect  with positive definition '
#ifdef DFS
     call print_maxmin_six(qp  ,ib*jbw          ,ntr*levs_  ,1,ntr*levs_ ,'qp   input')
     call print_maxmin_six(ut  ,ib*jbw          ,levs   ,1,levs  ,'ut   input')
     call print_maxmin_six(pdot,ib*jbw          ,levs+1 ,1,levs+1,'pdot input')
#else
     call print_maxmin_six(qp  ,LONF2S*ntr*levs_    ,LATG2S ,1,LATG2S,'qp   input')
     call print_maxmin_six(ut  ,LONF2S*levs_    ,LATG2S ,1,LATG2S,'ut   input')
     call print_maxmin_six(pdot,LONF2S*(levs_+1),LATG2S ,1,LATG2S,'pdot input')
#endif
   endif
#endif
!
! first mass conserving interpolation from reduced grid to full grid
!
   do jj = 1,jjend
     j1=2*jj-1     ! N.H.
     j2=2*jj       ! S.H.
#ifdef REDUCE_GRID
#ifdef MP
     lat=latdef(latstr(mype)+jj-1)
#else
     lat=latdef(jj)
#endif
     lonsd=lonfd(lat)
#else
     lat=latdef(jj)
     lonsd=nx
#endif /* REDUCE_GRID end */
     ! u,v at time step n (convert dynamics to SL grid)
     do k = 1,LEVSS
       do i = 1,lonsd
#ifdef DFS
         uulon(i,k,j1) = UU(i,jj      ,k) * (rrerth_*rbs2(jj))
         uulon(i,k,j2) = UU(i,jbw+1-jj,k) * (rrerth_*rbs2(jj))
         vvlon(i,k,j1) = VV(i,jj      ,k) * (rrerth_*sqrt(rbs2(jj)))
         vvlon(i,k,j2) = VV(i,jbw+1-jj,k) * (rrerth_*sqrt(rbs2(jj)))
#else
         uulon(i,k,j1) = UU(i      ,k,jj) * (rrerth_*rbs2(jj))
         uulon(i,k,j2) = UU(lonsd+i,k,jj) * (rrerth_*rbs2(jj))
         vvlon(i,k,j1) = VV(i      ,k,jj) * (rrerth_*sqrt(rbs2(jj)))
         vvlon(i,k,j2) = VV(lonsd+i,k,jj) * (rrerth_*sqrt(rbs2(jj)))
#endif
       enddo
     enddo
     ! at time step n-1 (convert dynamics to SL grid)
     do k = 1,ntr*LEVSS
       do i = 1,lonsd
#ifdef DFS
         qqlon(i,k,j1) = QQ(i,jj      ,k)
         qqlon(i,k,j2) = QQ(i,jbw+1-jj,k)
#else
         qqlon(i,k,j1) = QQ(i      ,k,jj)
         qqlon(i,k,j2) = QQ(lonsd+i,k,jj)
#endif
       enddo
     enddo
#ifdef SLDBG
     if( iope ) then
       print *,'jj,lat,lonsd in nislq_advect',jj,lat,lonsd
#ifdef DFS
       if(lat.eq.1) then
#else
       if(lat.eq.7) then
#endif
         call print_maxmin_six(uulon(1,1,j1),nx,LEVSS ,1,LEVSS ,               &
                                                         'u in dlamda/dt')
         call print_maxmin_six(vvlon(1,1,j1),nx,LEVSS ,1,LEVSS ,               &
                                                         'v in dphi/dt  ')
         call print_maxmin_six(qqlon(1,1,j1),nx,ntr*LEVSS,1,ntr*LEVSS,               &
                                                         'q in dq/dt    ')
       endif
     endif
#endif
#undef UU
#undef VV
#undef QQ
#ifdef REDUCE_GRID
!
! reduced grid to full grid
!
     call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,uulon(1,1,j1))
     call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,vvlon(1,1,j1))
     call cyclic_cell_intpx(ntr*LEVSS,lonsd,lonfull,qqlon(1,1,j1))
!
     call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,uulon(1,1,j2))
     call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,vvlon(1,1,j2))
     call cyclic_cell_intpx(ntr*LEVSS,lonsd,lonfull,qqlon(1,1,j2))
#endif
!
     do k = 1,ntr*LEVSS
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
#ifdef SLDBG
   if( iope ) then
     print *,' finish cyclic_cell_massadvx'
     call print_maxmin_six(rrlon,lonfull*ntr*LEVSS,latpart,1,latpart,             &
                                                          'r advx 1st     ')
     call print_maxmin_six(vvlon,lonfull*LEVSS ,latpart,1,latpart,             &
                                                          'vv before we2ns')
   endif
#endif
! ---------------------------------------------------------------------
! mpi para from horizontal full grid to meridional full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!
   call nislq_transpose_we2ns(vvlon,vvlat,LEVSS ,nsize)
   call nislq_transpose_we2ns(qqlon,qqlat,ntr*LEVSS,nsize)
   call nislq_transpose_we2ns(rrlon,rrlat,ntr*LEVSS,nsize)
#ifdef SLDBG
   if( iope ) then
     print *,' nislq transport from we to ns '
     call print_maxmin_six(vvlat,latfull*LEVSS ,mylonlen,1,mylonlen,'v we2ns')
     call print_maxmin_six(qqlat,latfull*ntr*LEVSS,mylonlen,1,mylonlen,'q we2ns')
     call print_maxmin_six(rrlat,latfull*ntr*LEVSS,mylonlen,1,mylonlen,'r we2ns')
   endif
#endif
! ---------------------------------------------------------------------
! ------------------- in meridional great circle ----------------------
! ---------------------------------------------------------------------
#ifdef SLDBG
   if ( iope ) then
     print *,' nislq adv loop in y '
     print *,' mylonlen=',mylonlen
   endif
#endif
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
#ifdef SLDBG
   if( iope ) then
     call print_maxmin_six(rrlat,latfull*ntr*LEVSS,mylonlen,1,mylonlen,           &
                                                           'r advy 1st')
     call print_maxmin_six(qqlat,latfull*ntr*LEVSS,mylonlen,1,mylonlen,           &
                                                           'q advy 2nd')
   endif
#endif
! ----------------------------------------------------------------------
! mpi para from meridional direction to horizontal directory 
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon
!
   call nislq_transpose_ns2we(qqlat,qqlon,ntr*LEVSS,nsize)
   call nislq_transpose_ns2we(rrlat,rrlon,ntr*LEVSS,nsize)
#ifdef SLDBG
   if( iope ) then
     print *,' nislfv_advq transport from ns to we '
     call print_maxmin_six(qqlon,lonfull*ntr*LEVSS,latpart,1,latpart,'q ns2we 2nd')
     call print_maxmin_six(rrlon,lonfull*ntr*LEVSS,latpart,1,latpart,'r ns2we 1st')
   endif
#endif
! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!     print *,' nislq adv loop in x for last '
!
   do jj = 1,jjend
     j1=2*jj-1     ! N.H.
     j2=2*jj       ! S.H.
#ifdef REDUCE_GRID
#ifdef MP
     lat=latdef(latstr(mype)+jj-1)
#else
     lat=latdef(jj)
#endif
     lonsd=lonfd(lat)
#else
     lat=latdef(jj)
     lonsd=nx
#endif
!
! second set advection in x for the second of the pair
!
     call cyclic_cell_massadvx(LEVSS,ntr,deltim,                              &
                                          uulon(1,1,j1),qqlon(1,1,j1),mass)
     call cyclic_cell_massadvx(LEVSS,ntr,deltim,                              &
                                          uulon(1,1,j2),qqlon(1,1,j2),mass)
     do k = 1,ntr*LEVSS
       do i = 1,lonsd
         rrlon(i,k,j1) = 0.5 * ( qqlon(i,k,j1) + rrlon(i,k,j1) )
         rrlon(i,k,j2) = 0.5 * ( qqlon(i,k,j2) + rrlon(i,k,j2) )
       enddo
     enddo
#ifdef REDUCE_GRID
!
! full gird to reduced grid
!
     call cyclic_cell_intpx(ntr*LEVSS,lonfull,lonsd,rrlon(1,1,j1))
     call cyclic_cell_intpx(ntr*LEVSS,lonfull,lonsd,rrlon(1,1,j2))
#endif
!
! convert SL to dynamics grid
!
     do k = 1,ntr*LEVSS
       do i = 1,lonsd
#ifdef DFS
         qt(i,jj      ,k)=rrlon(i,k,j1)
         qt(i,jbw+1-jj,k)=rrlon(i,k,j2)
#else
         qt(i      ,k,jj)=rrlon(i,k,j1)
         qt(lonsd+i,k,jj)=rrlon(i,k,j2)
#endif
       enddo
     enddo
   enddo
#ifdef MP
!
! transpose x-full to z-full
!
#ifdef DFS
   call mpyz2xy(qt,iba,levsp,qtp,ib,levs,jbw,ntr)
#else
   call mpnk2nx(qt,lonf2_,ntr*LEVSS,qtp,lonf2p_,ntr*levs_,latg2p_,levsp_,levs_,       &
                                                                  1,1,ntr)
#endif
#define QT qtp
#else
#define QT qt
#endif /* MP end */
! --------------------------------------------------------------
! ----------- compute vertical advection and total ------------
! --------------------------------------------------------------
#ifdef DFS
   do j = 1,jbw
#else
   do j = 1,jjend
#endif
#ifdef REDUCE_GRID
#ifdef MP
     lonsd=lonfdp(j,mype)*2
#else
     lonsd=lonfd(latdef(j))*2
#endif
#else
#ifdef DFS
     lonsd=ib
#else
     lonsd=LONF2S
#endif
#endif /* REDUCE_GRID end */
!
!    pressure (top to bottom)
!
     do k = 1,levs_+1
       do i = 1,lonsd
#ifdef HYBRID
         ppi(i,k)=ak5(k)+bk5(k)*pt(i,j)
#else
#ifdef DCMIP
         ppi(i,k)=ptop*1.e-3+si(levs_+2-k)*(pt(i,j)-ptop*1.e-3)  ! cb
#else
         ppi(i,k)=si(levs_+2-k)*pt(i,j)
#endif /* DCMIP end */
#endif /* HYBRID end */
#ifdef DFS
         pdot2(i,k)=pdot(i,j,k)
#else
         pdot2(i,k)=pdot(i,k,j)
#endif
       enddo
     enddo
!
     do k = 1,ntr*levs_
       do i = 1,lonsd
#ifdef DFS
         qtn(i,k)=QT(i,j,k)
#else
         qtn(i,k)=QT(i,k,j)
#endif
       enddo
     enddo
#undef QT
#ifdef SLDBG
     if( iope ) then
       if(j.eq.1) call print_maxmin_seven(ppi,lonsd,LONF2S,lev+1,1,lev+1,'ppi')
     endif
#endif
!
! vertical advection with mass conserving positive advection
!
     call vertical_cell_advect(lonsd,LONF2S,lev,ntr,deltim,ppi,pdot2,qtn,mass)
!
!    q update at time step n+1 (bottom to top)
!
     do k = 1,levs_*ntr
       do i = 1,lonsd
#ifdef DFS
         q3(i,j,levs*ntr+1-k)=qtn(i,k)
#else
         q3(i,levs_*ntr+1-k,j)=qtn(i,k)
#endif
       enddo
     enddo
   enddo
!
#ifdef SLDBG
   if( iope ) then
     call print_maxmin_six(qqlon,lonfull*ntr*LEVSS,latpart,1,latpart,'q advx 2nd  ')
     call print_maxmin_six(rrlon,lonfull*ntr*LEVSS,latpart,1,latpart,'r after mean')
#ifdef DFS
     call print_maxmin_six(q1   ,ib*jbw       ,levs*ntr   ,1,levs*ntr   ,'q1 input    ')
     call print_maxmin_six(q3   ,ib*jbw       ,levs*ntr   ,1,levs*ntr   ,'q3 output   ')
#else
     call print_maxmin_six(q1   ,LONF2S*levs_*ntr ,LATG2S ,1,LATG2S ,'q1 input    ')
     call print_maxmin_six(q3   ,LONF2S*levs_*ntr ,LATG2S ,1,LATG2S ,'q3 output   ')
#endif
   endif
#endif
!
   return
#ifdef DFS
#undef LONF2S
#undef LATG2S
#undef LEVSS
#undef levs_
#undef levh_
#endif
#endif /* NISLQ_MONO end */
   end subroutine nislq_chem_advect

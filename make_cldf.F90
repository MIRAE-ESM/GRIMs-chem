!-------------------------------------------------------------------------------
   subroutine make_cldf
!-------------------------------------------------------------------------------
   use comfphys, only     : xlon, xlat
   use comgrad, only      : slmskr, rhcl
   use dao_mod, only      : CLDF, SPHU, T
   use pressure_mod, only : GET_PCENTER, GET_PEDGE
   use cmn_size_mod
!-------------------------------------------------------------------------------
!
   integer,parameter    ::  mcld=3,nbin=100
!
   real                 ::  cldary(IIPAR,LLPAR)
   real                 ::  prsi(IIPAR,LLPAR+1)
   real                 ::  prsl(IIPAR,LLPAR)
   real                 ::  rhcld(IIPAR,nbin,mcld)
   real                 ::  qlat(IIPAR,LLPAR)
   real                 ::  tlat(IIPAR,LLPAR)
!
   ! initialize
   cldary = 0d0
   prsi = 0d0
   prsl = 0d0
   rhcld = 0d0
   qlat = 0d0
   tlat = 0d0

   do lat = 1,LATLOOP

     ! pedge
     do lev = 1,LLPAR+1
     do lon = 1,LONLOOP
       prsi(lon,lev) = GET_PEDGE(lon,lat,lev) * 1d-1
     enddo
     enddo

     ! pcenter, q, t
     do lev = 1,LLPAR
     do lon = 1,LONLOOP
       prsl(lon,lev) = GET_PCENTER(lon,lat,lev) * 1d-1
       qlat(lon,lev) = SPHU(lon,lat,lev) * 1d-3
       tlat(lon,lev) = T(lon,lat,lev)
     enddo
     enddo

     ! make_rhcld
     call make_rhcld(IIPAR,slmskr(1,lat),xlon(1,lat),xlat(1,lat),rhcl,rhcld)

     ! make_cldf
     call make_cldary(IIPAR,LLPAR,qlat,tlat,prsi,prsl,cldary,xlat(1,lat),rhcld)

     ! CLDF
     do lev = 1,LLPAR
     do lon = 1,LONLOOP
       CLDF(lev,lon,lat) = cldary(lon,lev)
     enddo
     enddo

   enddo
!
   return
   end subroutine make_cldf
!
!-------------------------------------------------------------------------------
   subroutine make_rhcld(lons2,slmskr,rlon,rlat,rhcl,rhcld)
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S
   use constant, only  : pi_
!-------------------------------------------------------------------------------
!
   integer,parameter    ::  mcld=3,nseal=2,nbin=100,nlon=2,nlat=4
!
   real                 ::  rlon(LONF2S),rlat(LONF2S)
   real                 ::  slmskr(LONF2S)
   real,save            ::  xlim
!
!  array added for rh-cl calculation
!  indices for lon,lat,cld type(l,m,h), land/sea respectively
!  nlon=1-2, for eastern and western hemispheres
!  nlat=1-4, for 60n-30n,30n-equ,equ-30s,30s-60s
!  land/sea=1-2 for land(and seaice),sea
!
   real                 ::  rhcl (nbin,nlon,nlat,mcld,nseal)
   real                 ::  rhcla(nbin,nlon,mcld)
   real                 ::  rhcld(LONF2S,nbin,mcld)
   real                 ::  temp1, temp2
   real,save            ::  xlabdy(3),xlobdy(3)
!
!  xlabdy = lat bndry between tuning regions,+/- xlim for transition
!  xlobdy = lon bndry between tuning regions
!
   data xlabdy / 30.e0 , 0.e0 , -30.e0 /
   data xlobdy / 0.e0 , 180.e0 , 360.e0 /
   data xlim / 5.e0 /
!
! initialize local variables
!
   rhcla=0. ; rhcld=0.
!
!   generalized the computation and suitable for all grid ---
!                                             by h.-m.juang
   do i = 1,lons2
     isla = 1
     if (slmskr(i).lt.1.e0) isla = 2
     xlatpt = rlat(i) * 180.e0 / pi_
!
!  get rh-cld relation for this lat
!
     kla = 4
     do k = 1,3
       if (xlatpt.gt.xlabdy(k)) then
         kla = k
!
         exit
!
       end if
     enddo
!
     klap=0
     do k = 1,3
       xlnn = xlabdy(k)+xlim
       xlss = xlabdy(k)-xlim
       if (xlatpt.lt.xlnn.and.xlatpt.gt.xlss) then
         kla  = k
         klap = k+1
!
         exit
!
       endif
     enddo
!
     if( klap .eq. 0 ) then
       do kc = 1,mcld
         do lo = 1,nlon
           do nbi = 1,nbin
             rhcla(nbi,lo,kc) = rhcl(nbi,lo,kla,kc,isla)
           enddo
         enddo
       enddo
     else
!
!  linear transition between latitudinal regions...
!
       temp1=(xlatpt-xlss)/(xlnn-xlss)
       do kc = 1,mcld
         do lo = 1,nlon
           do nbi = 1,nbin
             temp2=(rhcl(nbi,lo,kla,kc,isla)-rhcl(nbi,lo,klap,kc,isla))
             rhcla(nbi,lo,kc) = temp1 * temp2 + rhcl(nbi,lo,klap,kc,isla)
           enddo
         enddo
       enddo
     endif
!
!  get rh-cld relation for this lon
!
     xlonpt = rlon(i) * 180.e0 / pi_
     lo=1
     if( xlonpt .gt. 180.0 ) lo=2
     ikn = 0
     do k = 1,3
       diflo = abs(xlonpt-xlobdy(k))
       if (diflo.lt.xlim) then
         ikn = k
         ilft = lo
         irgt = ilft + 1
         if (irgt.gt.nlon) irgt = 1
         xlft = xlobdy(ikn) - xlim
         xrgt = xlobdy(ikn) + xlim
!
         exit
!
       endif
     enddo
     if( ikn.eq.0 ) then
       do k = 1,mcld
         do nbi=1,nbin
           rhcld(i,nbi,k) = rhcla(nbi,lo,k)
         enddo
       enddo
     else
       do k = 1,mcld
         do nbi = 1,nbin
           rhcld(i,nbi,k) =                                                    &
           (rhcla(nbi,ilft,k)-rhcla(nbi,irgt,k))                               &
                  * (xlonpt-xrgt)/(xlft-xrgt)+rhcla(nbi,irgt,k)
         enddo
       enddo
     endif
!
   enddo
!
   return
   end subroutine make_rhcld
!
!-------------------------------------------------------------------------------
   subroutine make_cldary(imx,kmx,q,t,prsi,prsl,cldary,xlatrd,rhcld)
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use constant, only  : pi=>pi_,rd=>rd_,rv=>rv_,akapa_
   use comcd1
!-------------------------------------------------------------------------------
   integer,parameter    ::  mcld=3, nbin=100
   real,parameter       ::  eps=rd/rv, epsm1=rd/rv-1.0
   real                 ::  prsi(imx,kmx+1), prsl(imx,kmx)
   real                 ::  t(imx,kmx), q(imx,kmx)
   real                 ::  cldary(imx,kmx),xlatrd(imx)
   real                 ::  cr1(imx),cr2(imx)
!
!   rh-cld relationships for each point
!
   dimension rhcld(imx,nbin,mcld)
!
!   ptopc(k,l): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!       l=1,2 are low-lat (<45 degree) and pole regions)
!
!  workspace
!
   logical              ::  bitx(ILOTS),bit1
   real                 ::  rhrh (ILOTS,levs_)
   real                 ::  prsly(ILOTS,levs_)
   real                 ::  dthdp(ILOTS,levs_)
   real                 ::  theta(ILOTS,levs_)
   real                 ::  ptop1(ILOTS,4)
   integer              ::  kcut (ILOTS),kbase(ILOTS)
   integer              ::  ksave(ILOTS)
!-------------------------------------------------------------------------------
!
! initialize local variables
!
   rhrh =0.;  prsly=0.;  dthdp=0.;  theta=0.
   ptop1=0.;  kcut =0.;  kbase=0.;  ksave=0.
!
!   begin here 
!
   kdim=kmx
   kdimp=kmx+1
   levm1=kmx-1
   levm2=kmx-2
!
!  find top pressure for each cloud domain
!
   do k = 1,4
     do i = 1,imx
       fac = max (0.0e0, 4.0e0*abs(xlatrd(i))/pi-1.0e0)
       ptop1(i,k) = ptopc(k,1) + (ptopc(k,2)-ptopc(k,1)) * fac
     enddo
   enddo
!
!  low cloud top sigma level, computed for each lat cause
!       domain definition changes with latitude...
   klow=kdim
   do k = kdim,1,-1
     do i = 1,imx
       if (prsi(i,k)/prsi(i,1) .lt. ptop1(i,2)*1.0e-3) klow = min(klow,k)
     enddo
   enddo
!
!  potential temp and layer relative humidity
!  
   do k = 1,kdim
     do i = 1,imx 
       cldary(i,k) = 0.0e0
       prsly(i,k) = prsl(i,k) * 10.0e0
       exnr = (prsly(i,k)*0.001e0) ** (-akapa_)
       theta(i,k) = exnr * t(i,k)
       es = fpvs0(t(i,k))
       qs = eps * es / (prsl(i,k) + epsm1*es)
       rhrh(i,k) = max (0.0e0, min (1.0e0, q(i,k)/qs))
     enddo
   enddo
!
!   potential temp lapse rate
!  
   do k = 1,levm1
     do i = 1,imx
       dthdp(i,k) = (theta(i,k+1) - theta(i,k)) /(prsly(i,k+1) - prsly(i,k))
     enddo
   enddo
! 
!     find the stratosphere cut off layer for high cloud. it
!      is assumed to be above the layer with dthdp less than
!      -0.25 in the high cloud domain (from looking at 1 case).
!  
   do i = 1,imx
     kcut(i) = levm2
   enddo
   do k = klow+1,levm2
     bit1 = .false.
     do i = 1,imx
       if (kcut(i).eq.levm2 .and. prsly(i,k).le.ptop1(i,3) .and.               &
               dthdp(i,k).lt.-0.25e0) then
         kcut(i) = k
       end if  
       bit1    = bit1 .or. kcut(i).eq.levm2
     enddo
     if (.not. bit1) exit
   enddo
!
   do klev = 1,levm2
     do i = 1,imx
       kbase(i)=0
       bitx(i)=.false.
     enddo
     do kc = mcld,1,-1
       do i = 1,imx 
         if(prsly(i,klev).ge.ptop1(i,kc+1)) kbase(i)=kc
       enddo
     enddo
     nx=0
     nhalf=(nbin+1)/2
     do i = 1,imx
       if(kbase(i).le.0.or.klev.gt.kcut(i)) then
         cldary(i,klev)=0.
       elseif(rhrh(i,klev).le.rhcld(i,1,kbase(i))) then
         cldary(i,klev)=0.
       elseif(rhrh(i,klev).ge.rhcld(i,nbin,kbase(i))) then
         cldary(i,klev)=1.
       else
         bitx(i)=.true.
         ksave(i)=nhalf
         nx=nx+1
       endif
     enddo 
     do while(nx.gt.0)
       nhalf=(nhalf+1)/2
       do i = 1,imx
         if(bitx(i)) then
           crk=rhrh(i,klev)
           cr1(i)=rhcld(i,ksave(i),kbase(i))
           cr2(i)=rhcld(i,ksave(i)+1,kbase(i))
           if(crk.le.cr1(i)) then
             ksave(i)=max(ksave(i)-nhalf,1)
           elseif(crk.gt.cr2(i)) then
             ksave(i)=min(ksave(i)+nhalf,nbin-1)
           else
             cldary(i,klev)=0.01*(ksave(i)+(crk-cr1(i))/(cr2(i)-cr1(i)))
             bitx(i)=.false.
             nx=nx-1
           endif
         endif
       enddo
     enddo
   enddo
!
   return
   end subroutine make_cldary

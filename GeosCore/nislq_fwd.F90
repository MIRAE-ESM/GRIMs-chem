!-------------------------------------------------------------------------------
   subroutine cyclic_cell_massadvx_fwd(levs,ncld,delt,uc,qq,mass)
!-------------------------------------------------------------------------------
!
! compute local positive advection with mass conservation
! qq is advected by uc which is in radiance/sec to next position
!
! author: hann-ming henry juang 2008
!
   use nislq, only : lonfull, ggloni, cyclic_cell_ppm_intp
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real   , parameter                    ::  fa1 = 9./16.                     ,&
                                             fa2 = 1./16.
   integer                               ::  levs,ncld,mass
   integer                               ::  i,k,im,n
   real                                  ::  delt, rm2, dist, sc
   real   , dimension(lonfull,levs)      ::  uc
   real   , dimension(lonfull,levs,ncld) ::  qq
   real   , dimension(lonfull)           ::  next,da,dxfact
   real   , dimension(lonfull+1)         ::  xnext,uint
!
! preparations ---------------------------
!
! x is equal grid spacing, so location can be specified by grid point number
!
   im = lonfull
!
   do k = 1,levs
!
! 4th order interpolation from mid point to cell interfaces
!
     do i = 3,im-1
       uint(i)=fa1*(uc(i,k)+uc(i-1,k))-fa2*(uc(i+1,k)+uc(i-2,k))
     enddo
     uint(2)=fa1*(uc(2,k)+uc(1 ,k))-fa2*(uc(3,k)+uc(im  ,k))
     uint(1)=fa1*(uc(1,k)+uc(im,k))-fa2*(uc(2,k)+uc(im-1,k))
     uint(im+1)=uint(1)
     uint(im  )=fa1*(uc(im,k)+uc(im-1,k)) -fa2*(uc(1,k)+uc(im-2,k))
!
! compute next positions of cell interfaces
!
     do i = 1,im+1
       dist     = uint(i) * delt
       xnext(i) = ggloni(i) + dist
     enddo
!      
     if( mass.eq.1 ) then
       do i = 1,im
         dxfact(i) = (ggloni(i+1)-ggloni(i)) / (xnext(i+1)-xnext(i))
       enddo
     endif
!
!  mass positive advection
!
     sc=ggloni(im+1)-ggloni(1)
     do n = 1,ncld
       da(1:im) = qq(1:im,k,n)
       if(mass.eq.1) da(1:im) = da(1:im) * dxfact(1:im)
       call cyclic_cell_ppm_intp(xnext,da,ggloni,next,im,im,im,sc)
       qq(1:im,k,n) = next(1:im)
     enddo
!
   enddo
!
   return
   end subroutine cyclic_cell_massadvx_fwd
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_massadvy_fwd(levs,ncld,delt,vc,qq,mass)
!-------------------------------------------------------------------------------
!
! compute local positive advection with mass conserving
! qq will be advect by vc to next location with delt
!
! author: hann-ming henry juang 2007
!
   use nislq, only : lathalf, latfull, gglati, cyclic_cell_ppm_intp
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real   , parameter                    ::  fa1 = 9./16.                     ,&
                                             fa2 = 1./16.
   integer                               ::  levs,ncld,mass                   ,&
                                             n,k,j,jm,jm2
   real                                  ::  delt,sc
   real   , dimension(latfull,levs)      ::  vc
   real   , dimension(latfull,levs,ncld) ::  qq
   real   , dimension(latfull)           ::  var,da,next,dyfact
   real   , dimension(latfull+1)         ::  ynext,dist
!
! preparations ---------------------------
!
   jm   = lathalf
   jm2  = latfull
!
   do k = 1,levs
!
     do j = 1,jm
       var(j)      =-vc(j   ,k) * delt
       var(j+jm)   = vc(j+jm,k) * delt
     enddo
!
     do j = 3,jm2-1
       dist(j)=fa1*(var(j)+var(j-1))-fa2*(var(j+1)+var(j-2))
     enddo
     dist(2)=fa1*(var(2)+var(1  ))-fa2*(var(3)+var(jm2  ))
     dist(1)=fa1*(var(1)+var(jm2))-fa2*(var(2)+var(jm2-1))
     dist(jm2+1)=dist(1)
     dist(jm2  )=fa1*(var(jm2)+var(jm2-1))-fa2*(var(1)+var(jm2-2))
!
     do j = 1,jm2+1
       ynext(j) = gglati(j) + dist(j)
     enddo
!
     if( mass.eq.1 ) then
       do j = 1,jm2
         dyfact(j) = (gglati(j+1)-gglati(j)) / (ynext(j+1)-ynext(j))
       enddo
     endif
!
! advection all in y
!
     sc=gglati(jm2+1)-gglati(1)
     do n = 1,ncld
       da(1:jm2) = qq(1:jm2,k,n)
       if( mass.eq.1 ) da(1:jm2) = da(1:jm2) * dyfact(1:jm2)
       call cyclic_cell_ppm_intp(ynext,da,gglati,next,jm2,jm2,jm2,sc)
       qq(1:jm2,k,n) = next(1:jm2)
     enddo
! 
   enddo
!
   return
   end subroutine cyclic_cell_massadvy_fwd
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vertical_cell_advect_fwd(lons,londim,levs,ncld,deltim,           &
                                       ppi,wwi,qql,mass)
!
   use nislq, only : vertical_cell_ppm_intp
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                               ::  lons,londim,levs,ncld,i,k,n,mass
   real                                  ::  deltim
   real   , dimension(londim,levs+1)     ::  ppi,wwi
   real   , dimension(londim,levs,ncld)  ::  qql
   real   , dimension(levs)              ::  dsfact
   real   , dimension(levs+1)            ::  ppii,ppia
   real   , dimension(levs,ncld)         ::  rqnn,rqda
!
   do i = 1,lons
     do k = 1,levs+1
       ppii(k)=ppi(i,k)
       ppia(k)=ppi(i,k)+wwi(i,k)*deltim
     enddo
!
     if( mass.eq.1) then
       do k = 1,levs
         dsfact(k)=(ppii(k)-ppii(k+1))/(ppia(k)-ppia(k+1))
       enddo
     endif
!
     do n = 1,ncld                              !hmhj nisl
       do k = 1,levs
         rqda(k,n) = qql(i,k,n)
       enddo
     enddo
     if( mass.eq.1 ) then
       do n = 1,ncld
         do k = 1,levs
           rqda(k,n) = rqda(k,n) * dsfact(k)
         enddo
       enddo
     endif
     call vertical_cell_ppm_intp(ppia,rqda,ppii,rqnn,levs,ncld)
     do n = 1,ncld
       do k = 1,levs
         qql(i,k,n)=rqnn(k,n)
       enddo
     enddo
   enddo
!
   return
   end subroutine vertical_cell_advect_fwd
!-------------------------------------------------------------------------------

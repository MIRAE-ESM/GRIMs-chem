#include <define.h>
   subroutine cnv_gaussian_lat(lats,gaulat)   
!-------------------------------------------------------------------------------
!
! - This file contains "sph_legendre_polinomial"
! - implicit double precision (a-h,o-z)                                       
!
   real ::  gaulat(lats)                                                    
!-------------------------------------------------------------------------------
   lat2 = lats / 2                                                           
   eps = 1.0e-14                                                             
   si = 1.0e+00                                                              
   k2 = 2 * lat2                                                             
   k1 = k2 - 1                                                               
   pi = atan(si) * 4.0e+00                                                   
!   dradz = pi / 360.0e+00                                                    
   dradz = pi / (4.*lat2)
   rad = 0.0e+00                                                             
   do k=1,lat2                                                           
      drad = dradz                                                            
      !
      10  call sph_legendre_polinomial(k2,rad,p2)                                                    
      !
      20  p1 = p2                                                                 
      rad = rad + drad                                                        
      !
      call sph_legendre_polinomial(k2,rad,p2)                                                    
      !
      if (sign(si,p1) .eq. sign(si,p2)) goto 20                               
      if (drad .lt. eps) goto 30                                              
      rad = rad - drad                                                        
      drad = drad * 0.1e+00                                                  
      go to 10                                                                
      30  continue                                                                
      gaulat(k) = 0.5e+00 * pi - rad                                          
      gaulat(lats-k+1) = -gaulat(k)                                           
   enddo
   if (mod(lats,2) .eq. 1) gaulat(lat2+1) = 0.0                              
!
   return                                                                    
   end subroutine cnv_gaussian_lat
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine sph_legendre_polinomial(n,rad,p)
!-------------------------------------------------------------------------------
!  
! abstract: evaluates the unnormalized legendre polynomial                      
!   of specified degree at a given colatitude using a standard                  
!   recursion formula.  real arithmetic is used.                                
!                                                                               
! program history log:                                                          
!   1988-04-01  joseph sela                                                       
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call sph_legendre_polinomial (n, rad, p)                                               
!   input argument list:                                                        
!     n        - degree of legendre polynomial.                                 
!     rad      - real colatitude in radians.                                    
!                                                                               
!   output argument list:                                                       
!     p        - real value of legendre polynomial.                             
!                                                                               
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in)   ::  n
   real,    intent(in)   ::  rad
   real,    intent(out)  ::  p
   real     ::  x,y1,y2,y3,g
   integer  ::  i
!-------------------------------------------------------------------------------
   x = cos(rad)
   y1 = 1.0
   y2=x
   do i=2,n
      g=x*y2
      y3=g-y1+g-(g-y1)/float(i)
      y1=y2
      y2=y3
   enddo
   p=y3
!
   return
   end subroutine sph_legendre_polinomial
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

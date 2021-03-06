!----------------------------------------------------
! Ellipse r(theta) in polar coordinate system
!----------------------------------------------------

      REAL(MK) FUNCTION polar_ellipse_r(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        polar_ellipse_r = SQRT(2.0_MK) * a * b / &
             SQRT( ( b**2-a**2 )* &
             COS(2.0_MK*(theta-phi)) + &
             a**2 + b**2 )
        
      END FUNCTION polar_ellipse_r
      
      
      REAL(MK) FUNCTION polar_ellipse_dr(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        polar_ellipse_dr = &
             ( polar_ellipse_r(a,b,theta,phi) ) ** 3.0 * &
             ( b**2-a**2 ) * SIN(2.0_MK*(theta-phi)) / &
             (2.0_MK*a**2*b**2)
        
        
      END FUNCTION polar_ellipse_dr
      
      
      REAL(MK) FUNCTION polar_ellipse_ddr(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        REAL(MK)                        :: r
        
        r =  polar_ellipse_r(a,b,theta,phi) 
        
        polar_ellipse_ddr = 3.0_MK * &
             ( polar_ellipse_dr(a,b,theta,phi) ) ** 2 / &
             r + r**3 * ( b**2-a**2 ) * COS(2.0_MK*(theta-phi)) / &
             a**2/b**2
        
        
      END FUNCTION polar_ellipse_ddr
      
     

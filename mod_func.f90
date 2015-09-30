module mod_func

	implicit none
	real(8) :: rampf
	real(8),parameter :: g = 9.87

contains

    function eti(x,y)

        use wave
        implicit  none
 
        real(8),intent(in):: x,y
        real(8) :: eti

        eti=amp*dcos(wk*(x*dcos(beta)+y*dsin(beta)) -w1*timerk)
        eti=rampf*eti

    end function eti


     function poxy(x,y,z)
        use wave

        implicit  none
        
        real*8,intent(in):: x,y,z
        real*8  poxy
        real*8  wkx,dum

        if(h.gt.0.0d0) then
            dum=amp*g/w1*dcosh(wk*(z+h))/dcosh(wk*h)
        else
            dum=amp*g/w1*dexp(wk*z)
        endif

        wkx=wk*( x*dcos(beta)+y*dsin(beta))
        poxy=dum*dsin(wkx-w1*timerk)
            

        !c             poxy = dcos(x)*dsin(y)*dexp(dsqrt(2.0d0)*z)*dcos(-w1*timerk)
        !c

        ! h   : water depth
        ! amp : amplitude of incident waves
        ! beta: incident angle
        ! w1  : angular frequency of incident waves
        ! v   : wave number in deep water
        ! wk  : wave number
        ! tper: wave period
        ! tstep: time step
        ! time : simulation time at each time step
        ! timerk: simulation time at each rk step
        ! rampf: ramp function for incident potential
        ! rampv: ramp function for damping
        !	   poxy=rampf*poxy

        end function poxy
!
! ----------------------------------------------

!ctime derivatives of incident wave potential
!c  

    function dpot(x,y,z)
        use wave
        implicit  none

        real*8,intent(in):: x,y,z
        real*8  dpot
        real*8  wkx,dum

        if(h.gt.0.0d0) then
            dum=-amp*g*dcosh(wk*(z+h))/dcosh(wk*h)
        else
            dum=-amp*g*dexp(wk*z)
        endif

        wkx=wk*( x*dcos(beta)+y*dsin(beta))
        dpot=dum*dcos(wkx-w1*timerk)
        !c
        !	   dpot=rampf*dpot
        !c        dpot = dcos(x)*dsin(y)*dexp(dsqrt(2.0d0)*z)*dcos(-w1*timerk)*-w1

        end function dpot


!

!C-------------------------------------------
!CSpacial Derivatives of incident wave potential
!CIORDER=0: for a current
!CIORDER=1: for the first order potential
!CIORDER=2: for the second order potential 

    subroutine  dinp(x,y,z,dpox,dpoy,dpoz)

        use wave
        implicit    none

        real*8,intent(in)::   x,y,z
        real*8,intent(out)::  dpox,dpoy,dpoz

        real*8 dum,wkx

        dum=amp*g/w1
        wkx=wk*(x*dcos(beta)+y*dsin(beta))

        if (h.gt.0.0d0) then           
        dpox= dum*wk*dcos(beta)*&
            &dcosh(wk*(z+h))/dcosh(wk*h)*dcos(wkx-w1*timerk)
        dpoy= dum*wk*dsin(beta)*&
            &dcosh(wk*(z+h))/dcosh(wk*h)*dcos(wkx-w1*timerk)
        dpoz= dum*wk*&
            &dsinh(wk*(z+h))/dcosh(wk*h)*dsin(wkx-w1*timerk)

        else

        dpox= dum*wk*dcos(beta)*dexp(wk*z)*dcos(wkx-w1*timerk)
        dpoy= dum*wk*dsin(beta)*dexp(wk*z)*dcos(wkx-w1*timerk)
        dpoz= dum*wk*dexp(wk*z)*dsin(wkx-w1*timerk)

        endif
!
!	   DPOX=RAMPF*DPOX
!	   DPOY=RAMPF*DPOY
!	   DPOZ=RAMPF*DPOZ

!C            DPOX = -DSIN(X)*DSIN(Y)
!C            DPOX = DPOX*DEXP(DSQRT(2.0D0)*Z)*DCOS(-W1*TimeRK)
!C            DPOY = DCOS(X)*DCOS(Y)*DEXP(DSQRT(2.0D0)*Z)*DCOS(-W1*TimeRK)

!C            DPOZ = DSQRT(2.0D0)
!C            DPOZ = DPOZ * DCOS(X)*DSIN(Y)*DEXP(DSQRT(2.0D0)*Z)
!C            DPOZ = DPOZ*DCOS(-W1*TimeRK)
!C      DPOX = 0
!C      DPOY = 0
!C      DPOZ = 0
!
        end subroutine  dinp 
        
        
!
! ======================================================
! Spacial potential and derivatives of fiction potential
!  Kind=0: for constant potential phi=1.0 
!  Kind=1: for a uniform current in the x-direction phi=x 
!  Kind=2: for a uniform current in the y-direction phi=y 
!  Kind=3: for a uniform current in the z-direction phi=z 
!
!
        SUBROUTINE  DINP0(KIND,X,Y,Z,PHI,DPOX,DPOY,DPOZ)
	  IMPLICIT    NONE
	  
	  INTEGER KIND
	  REAL*8,INTENT(IN)::   X,Y,Z
        REAL*8,INTENT(OUT)::  PHI,DPOX,DPOY,DPOZ

		  
		IF (KIND==0) THEN
		  PHI=1.0d0           
          DPOX=0.0d0
          DPOY=0.0d0
          DPOZ=0.0d0
        ELSE IF (KIND==1) THEN
		  PHI=X       
          DPOX=1.0d0
          DPOY=0.0d0
          DPOZ=0.0d0
        ELSE IF (KIND==2) THEN
		  PHI=Y           
          DPOX=0.0d0
          DPOY=1.0d0
          DPOZ=0.0d0
        ELSE IF (KIND==3) THEN
		  PHI=Z           
          DPOX=0.0d0
          DPOY=0.0d0
          DPOZ=1.0d0          
        ENDIF
!
        RETURN
        END SUBROUTINE  DINP0         
        

!C*********************************************************
!C*                                                       *
!C* CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
!C* FOR 6 NODE TRIANGULAR ELEMENTS                        *
!C*                                                       *
!C*********************************************************     

        SUBROUTINE SPFUNC6(SI,ETA,SF,DSF)
      	IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8)

        SF(1)=(1.0D0-SI-ETA)*(1.0D0-2.0D0*SI-2.0D0*ETA)
        SF(2)=SI *(2.0D0*SI -1.0D0)
        SF(3)=ETA*(2.0D0*ETA-1.0D0)
        SF(4)=4.0D0*SI*(1.0D0-SI-ETA)
        SF(5)=4.0D0*SI*ETA
        SF(6)=4.0D0*ETA*(1.0D0-SI-ETA)

        DSF(1,1)=4.0D0*SI+4.0*ETA-3.0D0
        DSF(1,2)=4.0D0*SI-1.0D0
        DSF(1,3)=0.0D0
        DSF(1,4)=4.0D0-8.0D0*SI-4.0D0*ETA
        DSF(1,5)=4.0D0*ETA
        DSF(1,6)=-4.0*ETA

        DSF(2,1)=4.0D0*SI+4.0D0*ETA-3.0D0
        DSF(2,2)=0.0D0
        DSF(2,3)=4.0D0*ETA-1.0D0
        DSF(2,4)=-4.0*SI
        DSF(2,5)=4.0*SI
        DSF(2,6)=4.0-8.0D0*ETA-4.0*SI
 

	RETURN
	END SUBROUTINE SPFUNC6




!C*********************************************************
!C*                                                       *
!C* CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
!C* FOR 8 NODE QUADRILATERIAL ELEMENTS                    *
!C*                                                       *
!C*********************************************************     

        SUBROUTINE SPFUNC8(SI,ETA,SF,DSF)
        IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8)

        SF(1)=-0.25D0*(1.0D0-SI)*(1.0D0-ETA)*(1.0D0+SI+ETA)
        SF(2)= 0.5D0 *(1.0D0-SI*SI)*(1.0-ETA)
        SF(3)= 0.25D0*(1.0D0+SI)*(1.0D0-ETA)*(SI-ETA-1.0D0)
        SF(4)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0D0+SI)
        SF(5)= 0.25D0*(1.0D0+SI)*(1.0D0+ETA)*(SI+ETA-1.0D0)
        SF(6)= 0.5D0 *(1.0-SI*SI)*(1.0D0+ETA)
        SF(7)=-0.25D0*(1.0D0-SI)*(1.0D0+ETA)*(SI-ETA+1.0D0)
        SF(8)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0-SI)            
      
        DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        DSF(1,2)=-SI*(1.0D0-ETA)     	
        DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        DSF(1,6)=-SI*(1.0D0+ETA)     
        DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        DSF(1,8)=-0.5D0*(1.0-ETA*ETA)

        DSF(2,1)= 0.25D0*(SI+2.0D0*ETA)*(1.0D0-SI)
        DSF(2,2)=-0.5D0*(1.0-SI*SI)
        DSF(2,3)= 0.25D0*(1.0D0+SI)*(2.0D0*ETA-SI)
        DSF(2,4)=-(1.0D0+SI)*ETA
        DSF(2,5)= 0.25*(1.0D0+SI)*(SI+2.0D0*ETA)
        DSF(2,6)= 0.5D0*(1.0D0-SI*SI)
        DSF(2,7)=-0.25D0*(1.0D0-SI)*(SI-2.0D0*ETA)
        DSF(2,8)=-(1.0D0-SI)*ETA   

	RETURN
	END SUBROUTINE SPFUNC8



!C*********************************************************
!C*                                                       *
!C* CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
!C* FOR 6 NODE TRIANGULAR ELEMENTS                        *
!C*                                                       *
!C*********************************************************     

        SUBROUTINE SPFUNC6_1(SI,ETA,SF,DSF,DDSF)
      	IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8),DDSF(3,8)

        SF(1)=(1.0D0-SI-ETA)*(1.0D0-2.0D0*SI-2.0D0*ETA)
        SF(2)=SI *(2.0D0*SI -1.0D0)
        SF(3)=ETA*(2.0D0*ETA-1.0D0)
        SF(4)=4.0D0*SI*(1.0D0-SI-ETA)
        SF(5)=4.0D0*SI*ETA
        SF(6)=4.0D0*ETA*(1.0D0-SI-ETA)
! ------------------------------------------
!
        DSF(1,1)=4.0D0*SI+4.0*ETA-3.0D0
        DSF(1,2)=4.0D0*SI-1.0D0
        DSF(1,3)=0.0D0
        DSF(1,4)=4.0D0-8.0D0*SI-4.0D0*ETA
        DSF(1,5)=4.0D0*ETA
        DSF(1,6)=-4.0*ETA

        DSF(2,1)=4.0D0*SI+4.0D0*ETA-3.0D0
        DSF(2,2)=0.0D0
        DSF(2,3)=4.0D0*ETA-1.0D0
        DSF(2,4)=-4.0*SI
        DSF(2,5)=4.0*SI
        DSF(2,6)=4.0-8.0D0*ETA-4.0*SI

!  ---------------------------
        DDSF(1,1)= 4.D0
        DDSF(1,2)= 4.0D0     	
        DDSF(1,3)= 0.0D0
        DDSF(1,4)=-8.0D0
        DDSF(1,5)= 0.0D0
        DDSF(1,6)= 0.0D0     

        DDSF(2,1)= 4.0D0
        DDSF(2,2)= 0.0D0
        DDSF(2,3)= 4.0D0
        DDSF(2,4)= 0.0D0
        DDSF(2,5)= 0.0D0
        DDSF(2,6)=-8.0D0

        DDSF(3,1)= 4.0D0
        DDSF(3,2)= 0.0D0   	
        DDSF(3,3)= 0.0D0
        DDSF(3,4)=-4.0d0
        DDSF(3,5)= 4.0D0
        DDSF(3,6)=-4.0d0     

	  RETURN
	 END SUBROUTINE SPFUNC6_1




        SUBROUTINE SPFUNC8_1(SI,ETA,SF,DSF,DDSF)
        IMPLICIT  NONE

        REAL*8,INTENT(IN) :: SI,ETA
        REAL*8,INTENT(OUT):: SF(8),DSF(2,8),DDSF(3,8)
!
        SF(1)=-0.25D0*(1.0D0-SI)*(1.0D0-ETA)*(1.0D0+SI+ETA)
        SF(2)= 0.5D0 *(1.0D0-SI*SI)*(1.0-ETA)
        SF(3)= 0.25D0*(1.0D0+SI)*(1.0D0-ETA)*(SI-ETA-1.0D0)
        SF(4)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0D0+SI)
        SF(5)= 0.25D0*(1.0D0+SI)*(1.0D0+ETA)*(SI+ETA-1.0D0)
        SF(6)= 0.5D0 *(1.0-SI*SI)*(1.0D0+ETA)
        SF(7)=-0.25D0*(1.0D0-SI)*(1.0D0+ETA)*(SI-ETA+1.0D0)
        SF(8)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0-SI)            
! ------------------------------------------------------------
      
        DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        DSF(1,2)=-SI*(1.0D0-ETA)     	
        DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        DSF(1,6)=-SI*(1.0D0+ETA)     
        DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        DSF(1,8)=-0.5D0*(1.0-ETA*ETA)
!
        DSF(2,1)= 0.25D0*(SI+2.0D0*ETA)*(1.0D0-SI)
        DSF(2,2)=-0.5D0*(1.0-SI*SI)
        DSF(2,3)= 0.25D0*(1.0D0+SI)*(2.0D0*ETA-SI)
        DSF(2,4)=-(1.0D0+SI)*ETA
        DSF(2,5)= 0.25*(1.0D0+SI)*(SI+2.0D0*ETA)
        DSF(2,6)= 0.5D0*(1.0D0-SI*SI)
        DSF(2,7)=-0.25D0*(1.0D0-SI)*(SI-2.0D0*ETA)
        DSF(2,8)=-(1.0D0-SI)*ETA   
! --------------------------------------------------------------

        DDSF(1,1)=0.5D0*(1.0D0-ETA)
        DDSF(1,2)=-(1.0D0-ETA)     	
        DDSF(1,3)= 0.5D0*(1.0D0-ETA)
        DDSF(1,4)= 0.0D0
        DDSF(1,5)= 0.5D0*(1.0D0+ETA)
        DDSF(1,6)=-(1.0D0+ETA)     
        DDSF(1,7)= 0.5D0*(1.0D0+ETA)
        DDSF(1,8)= 0.0D0

        DDSF(2,1)= 0.50D0*(1.0D0-SI)
        DDSF(2,2)= 0.0D0
        DDSF(2,3)= 0.50D0*(1.0D0+SI)
        DDSF(2,4)=-(1.0D0+SI)
        DDSF(2,5)= 0.50D0*(1.0D0+SI)
        DDSF(2,6)= 0.0D0
        DDSF(2,7)= 0.50D0*(1.0D0-SI)
        DDSF(2,8)=-(1.0D0-SI) 

        DDSF(3,1)= 0.25D0*(-2.0D0*SI-ETA*2.0D0+1.0D0)
        DDSF(3,2)= SI     	
        DDSF(3,3)= 0.25D0*(-2.0D0*SI+ETA*2.0D0-1.0D0)
        DDSF(3,4)= -ETA
        DDSF(3,5)= 0.25D0*(2.0D0*SI+ETA*2.0D0+1.0D0)
        DDSF(3,6)=-SI     
        DDSF(3,7)= 0.25D0*(2.0D0*SI-ETA*2.0D0-1.0D0)
        DDSF(3,8)= ETA

        RETURN
	  END SUBROUTINE SPFUNC8_1
       
end module



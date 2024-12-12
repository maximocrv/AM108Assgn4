      SUBROUTINE ELMT01(D,UL,XL,STRL,EPSL,QL,IX,TL,S,P,VEL,ACCEL,
     1 NDF,NDM,NST,NS,NQ,ISW)
C
C---- MOCK ELEMENT ROUTINE
C
c---- Use double precision variables, functions etc.
c---- Give following statement at beginning of each subroutine:
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      CHARACTER*8 HEAD
      CHARACTER*1 O
      COMMON /CHEAD/ O,HEAD(20)
      COMMON /CDATA/ NUMNP,NUMEL,NUMMAT,NEN,NEQ,IPR,NSDM,NQDM,NQUAD
      COMMON /ELDATA/ DM,NELMT,MA,MCT,IEL,NEL
      DIMENSION D(*),UL(NDF,*),XL(NDM,*),IX(*),TL(*),S(NST,*),P(*)
      DIMENSION STRL(NS,*),EPSL(NS,*),QL(NQ,*),VEL(NDF,*),ACCEL(NDF,*)
      DIMENSION B(4,18), BTRANS(18,4), DMAT(4,4), SHP(3,9), SG(16) 
      DIMENSION TG(16), WG(16)
C
C---- DESCRIPTION OF THE PARAMETER LIST AND SOME OF COMMON BLOCK VARIABLES
C
C---- D(*)             = MATERIAL PROPERTY ARRAY
C---- UL(NDF,NEL)      = NODAL DISPLACEMENTS
C---- VEL(NDF,NEL)     = NODAL VELOCITIES
C---- ACCEL(NDF,NEL)   = NODAL ACCELERATIONS
C---- XL(NDM,NEL)      = NODAL COORDINATES
C---- TL(NEL)          = NODAL TEMPERATURES
C---- IX(NEL)          = CONNECTIVITY ARRAY
C---- STRL(NS,NQUAD)   = ELEMENT STRESSES AT INTEGRATION AND OUTPUT POINTS
C---- EPSL(NS,NQUAD)   = ELEMENT STRAINS AT INTEGRATION AND OUTPUT POINTS
C---- QL(NQ,NQUAD)     = INTERNAL VARIABLES AT INTEGRATION AND OUTPUT POINTS
C---- S(NST,NST)       = TANGENT STIFFNESS MATRIX (ISW = 3) OR
C----                    CONSISTENT MASS MATRIX (ISW = 5)
C---- P(NST)           = ELEMENT FORCE ARRAY (ISW = 6) OR
C----                    LUMPED MASS MATRIX (ISW = 5)
C---- HEAD(20)         = HEADER FOR PRINTOUTS
C---- THE FOLL VARIABLES HAVE BEEN DEFINED ELSEWHERE IN THE PROGRAM
C---- NDF              = NUMBER OF DEGREES OF FREEDOM PER NODE
C---- NDM              = NUMBER OF SPATIAL DIMENSIONS
C---- NST              = DIMENSION OF STIFFNESS ARRAY (NDF*NEN)
C---- NS               = DIMENSION OF STRESS OR STRAIN ARRAY (NSDM)
C---- NQ               = DIMENSION OF INTERNAL VARIABLE ARRAY (NQDM)
C---- NQUAD            = NUMBER OF INTEGRATION AND OUTPUT QUADRATURE POINTS
C---- NELMT            = ELEMENT NUMBER CORRESPONDING TO PRESENT ELEMENT
C---- NEL              = NUMBER OF NODES IN THE PRESENT ELEMENT
C---- NEN              = MAXIMUM NUMBER OF NODES PER ELEMENT IN THE MESH
C---- MA               = MATERIAL SET NUMBER TO WHICH THE PRESENT 
C----                    ELEMENT BELONGS
C---- ISW              = SWITCH WHICH IS SET BY THE CALLING MACRO AND
C                        WHICH CONTROLS TO WHICH MODULE OF THE ELEMENT ROUTINE
C                        THE EXECUTION IS TO BE TRANSFERRED:
C                        ISW = 1 -> READ IN MATERIAL PROPERTIES
C                        ISW = 2 -> UNUSED ISW
C                        ISW = 3 -> COMPUTE ELEMENT STIFFNESS MATRIX
C                        ISW = 4 -> PRINTOUT ELEMENT VARIABLES (E.G. STRESSES)
C                        ISW = 5 -> COMPUTE MASS MATRIX
C                        ISW = 6 -> COMPUTE ELEMENT FORCE ARRAY
C                        ISW = 7 -> UNUSED ISW
C                        ISW = 8 -> UNUSED ISW
C                        ISW = 9 -> UPDATE ELEMENT VARIABLES
C
C---- OUTLINE OF THE SUBROUTINE
C
C---- GO TO CORRECT ARRAY PROCESSOR
      GO TO (1,2,3,4,5,6,7,8,9), ISW
C
1     CONTINUE
C
C---- INPUT MATERIAL PROPERTIES
C---- THIS MODULE IS ACCESSED WHEN MACRO 'MATE' IS EXECUTED
C
C---- TASKS TO BE PERFORMED IN THIS MODULE:
C
C---- 1) READ INTO ARRAY 'D' MATERIAL PROPERTIES FOR THE ELEMENT
c----    (Make sure to enter above data after MATE card following the
c----     format you have specified in the READ statement here)
C---- 2) PRINT OUT AN ECHO OF THE INPUT MATERIAL DATA

      ! Read the material properties from unit 5
      read(5,*) D(1), D(2), D(3), D(4), D(5), D(6)
      
      ! Print out an echo of the input material data
      WRITE(6,*) '--- Material Data Input Echo ---'
      WRITE(6,*) 'D(1): Body Force X or R component = ', D(1)
      WRITE(6,*) 'D(2): Body Force Y or Z component = ', D(2)
      WRITE(6,*) 'D(3): Mode Switch = ', D(3)
      WRITE(6,*) 'D(4): Mass Density (0.0 for static problems) = ', D(4)
      WRITE(6,*) 'D(5): Young''s Modulus = ', D(5)
      WRITE(6,*) 'D(6): Poisson''s Ratio = ', D(6)
      WRITE(6,*) '--------------------------------'
C
      RETURN
C
      
      
2     CONTINUE
C---- UNUSED ISW
      RETURN
C
      
      
      
3     CONTINUE
C
C---- COMPUTE ELEMENT STIFFNESS MATRIX
C---- THIS MODULE IS ACCESSED WHEN MACROS 'TANG' AND 'UTAN' ARE EXECUTED
C
C---- TASKS TO BE PERFORMED IN THIS MODULE:
C
C---- 1) COMPUTE THE ELEMENT STIFFNESS MATRIX
C---- 2) RETURN IT IN ARRAY 'S(NST,NST)'
C
      MODE = D(3)
     
      ! Integration order
      IF (NEL == 4) THEN 
      NQP = 2  !(2*2 INTEGRATION ORDER)
      ELSEIF (NEL > 4) THEN 
      NQP = 3  !(3*3 INTEGRATION ORDER) 
      END IF 
      
C---- THE ROUTINE PGAUS2 
C---- WILL RETURN THE NUMBER OF QUADRATURE POINTS 'LINT' (4 FOR NQP=2 AND
                                                          !9 FOR NQP=3)
C---- AND THE NATURAL COORDINATES 'SG(LINT)' AND 'TG(LINT)' AND THE WEIGHTS
C---- 'WG(LINT)' CORRESPONDING TO THE ABOVE 'LINT' GAUSS POINTS.
      CALL PGAUS2(NQP,LINT,SG,TG,WG)
      
      ! LOOP OVER THE INTEGRATION POINTS:
      DO 300 L=1,LINT
      
      ! AT EVERY INTEGRATION POINT, COMPUTE SHAPE FUNCTIONS AND
C---- DERIVATIVES 'SHP' AND JACOBIAN 'XSJ'
      CALL SHAP2D(SG(L),TG(L),XL,SHP,XSJ,NDM,NEL,IX,.TRUE.)
      
      ! SET UP THE STRAIN-DISPLACEMENT MATRIX 'B(NSDM,NST)'.
      CALL BMAT(B,SHP,XL,XJC,NEL,NDM,NST,MODE)

C---- SET UP THE MATERIAL MATRIX 'DMAT(NSDM,NSDM)' CORRESPONDING TO 
C---- PLANE STRAIN OR PLANE STRESS OR AXISYMM 
      CALL ELAS(D,DMAT,MODE)
      
C---- ADD CONTRIBUTION TO THE STIFFNESS MATRIX 'S(NST,NST)' ARISING
C---- FROM THE CURRENT GAUSS POINT L.
      
      do i = 1, NST           ! Loop over rows of S
          do j = 1, NST       ! Loop over columns of S
              do k = 1, 4     ! Loop over rows of B (transposed) and DMAT
                  do m = 1, 4 ! Loop over columns of DMAT and rows of B
                     S(i,j)=S(i,j)+XSJ*XJC*WG(L)*B(k,i)*DMAT(k,m)*B(m,j)
                  end do
              end do
          end do
      end do
      
300   continue 
      
      RETURN
      
      
      
      
4     CONTINUE
C
C---- OUTPUT ELEMENT VARIABLES SUCH AS STRESSES AND STRAINS
C---- THIS MODULE IS ACCESSED WHEN MACRO 'STRE' IS EXECUTED
C
C---- TASKS TO BE CARRIED OUT IN THIS MODULE:
C
C---- 1) PRINTOURT ELEMENT VARIABLES SUCH AS STRESSES, STRAINS AND
C----    INTERNAL VARIABLES (IF ANY) AT OUTPUT QUADRATURE POINTS
C
C---- SET UP SHAPE FUNCTIONS ONLY AT OUTPUT QUADRATURE POINT (NOTE LOGICAL
C---- FLAG IN following CALL STATEMENT IS FALSE.
      CALL SHAP2D(0.0,0.0,XL,SHP,XSJ,NDM,NEL,IX,.FALSE.)
      
      MODE = D(3)
      
      ! CALCULATE THE COORDINATES OF THE GP BY INTERPOLATING THE SHAPE FUNCTIONS
      XGAUS = 0.d0
      YGAUS = 0.d0
      
      DO 15 I=1,NEL

          XGAUS = XGAUS + SHP(3,I)*XL(1,I)
          YGAUS = YGAUS + SHP(3,I)*XL(2,I) 
          
15    CONTINUE 
      
C---- PRINT OUT ELEMENT NUMBER 'NELMT', LOCATION OF OUTPUT QUADRATURE
C---- POINT 'XGAUS' (OR 'RGAUS') AND 'YGAUS' (OR 'ZGAUS') AND VALUE OF 
C---- STRESSES AND STRAINS :
      WRITE(6,*) '--- Output Element Variables ---'
      WRITE(6,*) 'Element Number = ', NELMT
      
      IF (MODE == 1 .OR. MODE == 3) THEN
          WRITE(6,*) 'X Coordinate of the Gauss Point = ', XGAUS
          WRITE(6,*) 'Y Coordinate of the Gauss Point = ', YGAUS
      ELSE 
          WRITE(6,*) 'R Coordinate of the Gauss Point = ', XGAUS
          WRITE(6,*) 'Z Coordinate of the Gauss Point = ', YGAUS
      END IF 
      
      WRITE(6,*) 'Stress at Gauss Point = ', (STRL(I,NQUAD),I=1,NSDM)
      WRITE(6,*) 'Strain at Gauss Point = ', (EPSL(I,NQUAD),I=1,NSDM)
      WRITE(6,*) '--------------------------------'
      
C
      RETURN
C
5     CONTINUE ! WE DO NOT NEED TO DO THIS ONE
C
C---- COMPUTE THE CONSISTENT AND LUMPED MASS MATRICES
C---- THIS MODULE IS ACCESSESD WHEN MACROS 'LMAS' AND 'CMAS' ARE EXECUTED
C
C---- TASKS TO BE CARRIED OUT IN THIS MODULE:
C
C---- 1) COMPUTE CONSISTENT AND LUMPED ELEMENT MASS MATRICES
C---- 2) RETURN THEM IN ARRAYS 'S(NST,NST)' AND 'P(NST)', RESPECTIVELY
C
      RETURN
C
6     CONTINUE
C
C---- COMPUTE ELEMENT FORCE ARRAY
C---- THIS MODULE IS ACCESSED WHEN MACRO 'FORM' IS EXECUTED
C
C---- TASKS TO BE PERFORMED IN THIS MODULE:
C
C---- 1) COMPUTE NODAL FORCES ARISING FROM EXTERNALLY APPLIED DISTRIBUTED
C----    LOADS (GRAVITY LOADS, HYDROSTATIC PRESSURES, THERMAL LOADS, ETC.)
C---- 2) COMPUTE INTERNAL FORCES ARISING FROM (INITIAL) ELEMENT STRESS FIELD
C---- 3) SUBTRACT THE LATTER FROM THE FORMER. THIS IS THE ELEMENT FORCE
C----    VECTOR. IN NON-LINEAR PROBLEMS, THIS IS CALLED AS OUT-OF-BALANCE
C----    FORCES.
C---- 4) RETURN THE ELEMENT FORCES IN ARRAY 'P(NST)'
C
C---- REMARKS:
C
C---- 1) THE INTERNAL FORCES DUE TO ELEMENT STRESSES IS GIVEN BY
C----    INTEGRAL OF 'BTRANS*STRL' OVER THE VOLUME OF THE ELEMENT.
C
      MODE = D(3)

      ! Integration order
      IF (NEL == 4) THEN 
          NQP = 2  !(2*2 INTEGRATION ORDER)
      ELSEIF (NEL > 4) THEN 
          NQP = 3  !(3*3 INTEGRATION ORDER) 
      END IF 
      
      CALL PGAUS2(NQP,LINT,SG,TG,WG)
C---- NEXT LOOP OVER INTEGRATION POINTS:
      DO 600 L=1,LINT
C---- COMPUTE SHAPE FUNCTIONS AND DERIVATIVES AT CURRENT GAUSS POINTS
C---- (SEE ISW=3)
      CALL SHAP2D(SG(L),TG(L),XL,SHP,XSJ,NDM,NEL,IX,.TRUE.)
C---- WITH THIS SET UP STRAIN-DISPLACEMENT MATRIX 'B(NSDM,NST)'
      CALL BMAT(B,SHP,XL,XJC,NEL,NDM,NST,MODE)
C---- FIRST  ADD TO ELEMENT FORCE ARRAY 'P' CONTRIBUTION DUE TO
C---- BODY FORCES ARISING FROM CURRENT GAUSS POINT 'L'  
      
      DO I=1,NEL      
          P(2*I-1) = P(2*I-1) + XSJ*XJC*WG(L)*SHP(3,2*I-1)*D(1)
          P(2*I) = P(2*I) + XSJ*XJC*WG(L)*SHP(3,2*I)*D(2)
      END DO      
      
C---- FROM THE EXTERNAL FORCE SO OBTAINED SUBTRACT INTERNAL FORCES 
C---- DUE TO ELEMENT STRESSES AT CURRENT GAUSS POINT L           
  
      !P = P - XSJ*XJC*WG(L)*MATMUL(TRANSPOSE(B),STRL(:,L))
      
      ! Step: Perform the operation directly on P
      DO i = 1, 2 * NEL
          DO j = 1, NS
              P(i) = P(i) - XSJ * XJC * WG(L) * B(j,i) * STRL(j,L)
          END DO
      END DO
      
600   CONTINUE
C---- LOOP OVER GAUSS POINTS ENDS HERE
      RETURN
C
7     CONTINUE
C---- UNUSED ISW
      RETURN
C
8     CONTINUE
C---- UNUSED ISW
      RETURN
C
9     CONTINUE
C
C----  UPDATE ELEMENT VARIABLES AT INTEGRATION
C---- AND OUTPUT POINTS. THIS MODULE IS ACCESSED WHEN MACRO 'CEQS'
C---- IS EXECUTED.
C
C---- TASKS TO BE PERFORMED IN THIS MODULE:
C
C---- 1) UPDATE ELEMENT STRESSES, STRAINS AND INTERNAL VARIABLES
C---- 2) RETURN UPDATED QUANTITIES IN 'STRL', 'EPSL' AND 'QL', RESPECTIVELY.
C
C---- REMARKS:
C
C---- 1) FOR LINEAR PROBLEMS, STRESSES AND STRAINS ARE COMPUTED DIRECTLY
C----    FROM DISPLACEMENTS AND STORED IN ARRAYS 'STRL', 'EPSL'.
C---- 2) FOR NON-LINEAR PROBLEMS, STRESSES AND STRAINS CANNOT BE COMPUTED
C----    DIRECTLY FROM DISPLACEMENTS AND HAVE TO BE INTEGRATED ALONG
C----    THE ENTIRE LOADING PATH. IN THIS CASE, STRESSES, STRAINS AND ALSO
C----    INTERNAL VARIABLES NEED TO BE UPDATED AND STORED in 'STRL', 'EPSL',
C----    AND 'QL' RESPECTIVELY. 
C     

      MODE = D(3)
      LL = 0
      
      DO 900  KKK  = 1,2
C----  KKK = 1: Update Stresses at 3x3 G.P. (NEL>4) & 2x2 G.P. (NEL=4);  
c----  KKK = 2 : Update stresses at 1x1 GP
C                                             
C---- SET UP INTEGRATION POINTS AND WEIGHTS
C---- AS BEFORE, NQP=2  IF NEL =4
C----        AND NQP=3  IF NEL >4
          
      ! Integration order
      IF (NEL == 4) THEN 
          NQP = 2  !(2*2 INTEGRATION ORDER)
      ELSEIF (NEL > 4) THEN 
          NQP = 3  !(3*3 INTEGRATION ORDER) 
      END IF 
          
      IF (KKK == 1) THEN
          CALL PGAUS2(NQP,LINT,SG,TG,WG)
      ELSE
          CALL PGAUS2(1,LINT,SG,TG,WG)
      END IF
      
C---- NEXT LOOP OVER INTEGRATION POINTS
      DO 900 L=1,LINT
          
          ! initialize both EPSL and STRL
          
          
          
      LL = LL + 1
C---- SET UP SHAPE FUNCTIONS AND DERIVATIVES AT CURRENT GAUSS POINT
      CALL SHAP2D(SG(L),TG(L),XL,SHP,XSJ,NDM,NEL,IX,.TRUE.)
C---- THE SUBROUTINE SHAP2D WILL RETURN SHAPE FUNCTIONS AND DERIVATIVES
C---- AT CURRENT GAUSS POINT 'L'
C---- WITH THIS SET UP STRAIN DISPLACEMENT MATRIX 'B(NSDM,NST)'
      CALL BMAT(B,SHP,XL,XJC,NEL,NDM,NST,MODE)
C---- NEXT SET UP MATERIAL MATRIX 'DMAT(NSDM,NSDM)' CORRESPONDING TO 
C---- PLANE STRAIN OR PLANE STRESS OR AXISYMM
      CALL ELAS(D,DMAT,MODE)
C---- USING STRAIN DISPLACEMENT MATRIX 'B(NSDM,NST)' AND NODAL DISPLACEMENT
C---- ARRAY 'UL(NDF,NEL)', COMPUTE STRAINS AT CURRENT QUADRATURE POINT L
C---- AND STORE IT IN '(EPSL(I,LL),I=1,NSDM)'
C---- NOTE: LL = 1 to 4 correspond to 2x2 GP and  LL = 5 to 1x1 GP if NEL=4.
C---- NOTE: LL = 1 to 9 correspond to 3x3 GP and  LL = 10 to 1x1 GP if NEL>4.
      
      DO I = 1,NEL
          DO J = 1,4
              EPSL(J,LL)=EPSL(J,LL)+B(J,2*I-1)*UL(1,I)+B(J,2*I)*UL(2,I)
          END DO
      END DO
      
      IF (MODE == 1) THEN
          EPSL(1,LL) = 0.d0
      ELSEIF (MODE == 3) THEN
          EPSL(1,LL) = -D(6)/(1-D(6))*(EPSL(2,LL)+EPSL(3,LL))
      END IF
      
C---- NEXT USING 'DMAT(NSDM,NSDM)' AND THE STRAINS '(EPSL(I,LL),I=1,NSDM)'  
C---- COMPUTE STRESSES '(STRL(I,LL),I=1,NSDM)'
      STRL(:,LL) = STRL(:,LL)  + MATMUL(DMAT,EPSL(:,LL))        
      
900   CONTINUE
C---- END OF LOOP OVER INTEGRATION POINTS AND OUTPUT POINT
      
      RETURN
C
C---- FORMAT STATEMENTS:
C
      END


      ! SUBROUTINE TO COMPUTE THE STRAIN-DISPLACEMENT MATRIX
      SUBROUTINE BMAT(B,SHP,XL,XJC,NEL,NDM,NST,MODE)
      Implicit double precision (a-h,o-z)
C
C---- SET UP B-MATRIX FOR 2D PLANE STRAIN/PLANE STRESS/AXISYMM
      
      DIMENSION B(4,*),SHP(3,*),XL(NDM,*)
      PI = 3.141592653589793d0
      
C---- ALSO COMPUTE: XJC = 2*PI*RADIUS FOR AXISYMM (MODE=3)
C----                   = 1.0 FOR PLANE STRAIN/PLANE STRESS (MODE=1 OR 2) 
      IF (MODE == 1 .OR. MODE == 3) THEN
          XJC = 1.d0
      ELSEIF (MODE == 2) THEN
      ! Calculate RADIUS
          RADIUS = 0.d0
          DO I = 1, NEL
              RADIUS = RADIUS + SHP(3, I) * XL(1, I)
          END DO
          XJC = 2.d0 * PI * RADIUS
      END IF
      
       
C---- INITALIZE
      DO 10 I=1,NST
      DO 10 J=1,4
      B(J,I) = 0.d0
10    CONTINUE
C---- COMPLETE
c---- LEAVE FIRST ROW OF B (CORRESP TO EPS_33) AS ZEROS FOR PLANE STRAIN / PLANE STRESS
C---- FIRST ROW OF B SHOULD HAVE TERMS CORRESP TO EPS_THETA_THETA FOR AXISYMM CASE
C
      IF (MODE == 1 .OR. MODE == 3) THEN
      DO 11 I=1,NEL   
          B(2,2*I-1) = SHP(1,I)
          B(3,2*I) = SHP(2,I)
          B(4,2*I-1) = SHP(2,I)
          B(4,2*I) = SHP(1,I)
11    CONTINUE
      
      ELSEIF (MODE == 2) THEN
      DO 12 I=1,NEL   
          B(1,2*I-1) = SHP(3,I)/RADIUS
          B(2,2*I-1) = SHP(1,I)
          B(3,2*I) = SHP(2,I)
          B(4,2*I-1) = SHP(2,I)
          B(4,2*I) = SHP(1,I)
12    CONTINUE 
      END IF
      
      RETURN
      END
      
      
      ! SUBROUTINE TO CODE THE MATERIAL MATRIX
      SUBROUTINE ELAS(D,DMAT,MODE)
C
      Implicit double precision (a-h,o-z)
C---- SET UP D-MATRIX FOR 2D PLANE STRAIN/PLANE STRESS/AXISYMM
C
      DIMENSION D(*), DMAT(4,4)      
C
C---- INITALIZE
      DO 10 I=1,4
      DO 10 J=1,4
          DMAT(I,J) = 0.d0
10    CONTINUE                               
C---- COMPLETE
C---- LEAVE 1ST ROW  & 1st COLUMN OF DMAT (CORRESP TO S33) AS ZERO FOR PLANE STRESS
C
      IF (MODE == 1 .OR. MODE == 2) THEN 
          DMAT(1,1) = (1.d0-D(6))*D(5)/((1.d0+D(6))*(1.d0-2*D(6)))
          DMAT(1,2) = (D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          DMAT(1,3) = (D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          
          DMAT(2,1) = (D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          DMAT(2,2) = (1-D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          DMAT(2,3) = (D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          
          DMAT(3,1) = (D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          DMAT(3,2) = (D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          DMAT(3,3) = (1-D(6))*D(5)/((1+D(6))*(1-2*D(6)))
          
          DMAT(4,4) = ((1-2*D(6))/2)*(D(5)/((1+D(6))*(1-2*D(6))))
          
      ELSEIF (MODE == 3) THEN 
          DMAT(2,2) = D(5)/(1-D(6)**2)
          DMAT(2,3) = (D(6))*D(5)/(1-D(6)**2)
          
          DMAT(3,2) = (D(6))*D(5)/(1-D(6)**2)
          DMAT(3,3) = D(5)/(1-D(6)**2)
          
          DMAT(4,4) = ((1-D(6))/2)*(D(5)/(1-D(6)**2))
      END IF 
      
      RETURN
      END
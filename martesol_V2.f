C=======================================================================
      PROGRAM MARTESOL
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER DD,MM,YY,HR,MIN,K,FEOF,NP,KMAX,IFROT
      PARAMETER ( NP = 10000000 )
      REAL*8 PI,RADM,SEGS,LAT,COLAT(NP),LONG,ELONG(NP),H,Z(NP),JD,D,T,
     &       THETA,ZETA,XI,R(3),RME(3),RMM(3),RS(3),RSH(3),AC,PYTHAG,
     &       ALT,DMS,UAKM,RMME(3),RMMZ(3),RSLG(3),RSLT(3),DATENUM,EPSE,
     &       JDI,JDF,JD0,DT,LAM,HROT(3),OBL,OBL0,I,OM,FDT,THP,PROT,
     &       HROTEA(3),HROTI(3),HROTO(3),
     &       XPNE,YPNE,ZPNE,XPNI,YPNI,ZPNI,XPNEI,YPNEI,
     &       ZPNEI,XPNO,YPNO,ZPNO
      CHARACTER*200 ARG,ALGO,PATHOUT
      CHARACTER*10 CHDIA,CHHORA
      LOGICAL LROT1,LROT2
C-----------------------------------------------------------------------
C     PARAMETROS GENERALES
C-----------------------------------------------------------------------
      PI = 4.D0*DATAN(1.D0)
      RADM = 3396.19D0
      UAKM = 149597870.7D0
      EPSE = 23.43928D0*PI/180.D0
C-----------------------------------------------------------------------
C     LECTURA DE PARAMETROS EN LINEA DE COMANDOS
C-----------------------------------------------------------------------
      OPEN(UNIT=10,FILE='MARTESOLV2.IN',STATUS='UNKNOWN')
C-----------------------------------------------------------------------
C     INGRESO DE DATOS DE ENTRADA
C-----------------------------------------------------------------------
      READ(10,*) ALGO
      READ(10,*) ALGO
      READ(10,*) DD,MM,YY
      READ(10,*) ALGO
      READ(10,*) HR,MIN,SEGS
      CALL CAL2JD(DD,MM,YY,HR,MIN,SEGS,JDI)
      READ(10,*) ALGO
      READ(10,*) DD,MM,YY
      READ(10,*) ALGO
      READ(10,*) HR,MIN,SEGS
      CALL CAL2JD(DD,MM,YY,HR,MIN,SEGS,JDF)
      READ(10,*) ALGO
      READ(10,*) FDT
      READ(10,*) ALGO
      READ(10,*) IFROT
      READ(10,*) ALGO
      READ(10,*) OBL0
      READ(10,*) ALGO
      READ(10,*) ARG
      READ(10,*) ALGO
      READ(10,*) PATHOUT
C     ------------------------------------------------------------------
      LROT1 = IFROT.EQ.1
      LROT2 = IFROT.EQ.2
C-----------------------------------------------------------------------
      CALL CAL2JD(1,1,0,0,0,0.D0,JD0)
C-----------------------------------------------------------------------
C     ARCHIVO DE ENTRADA
C-----------------------------------------------------------------------
      OPEN(UNIT=11,FILE=TRIM(ARG))
C-----------------------------------------------------------------------
C     LECTURA DE DATOS EN ARCHIVO DE ENTRADA
C-----------------------------------------------------------------------
         K = 1
      FEOF = 0
      PRINT *, 'PROCESANDO ARCHIVO "'//TRIM(ARG)//'"'
      DO WHILE(FEOF.EQ.0)
C      -----------------------------------------------------------------
       READ(11,*,IOSTAT=FEOF) LONG,LAT,H
C      -----------------------------------------------------------------
       IF (FEOF.GT.0) THEN
        PRINT *,'REVISAR ARCHIVO DE ENTRADA'
        EXIT
       ELSE IF (FEOF.LT.0) THEN
        PRINT *,  'LISTO!'
        EXIT
       ELSE
C      -----------------------------------------------------------------
        IF (LAT.GE.0.D0) THEN
         COLAT(K) = 90.D0 - LAT
        ELSE
         COLAT(K) = 90.D0 + DABS(LAT)
        END IF
C       ----------------------------------------------------------------
        IF (LONG.GE.0.D0) THEN
         ELONG(K) = LONG
        ELSE
         ELONG(K) = 360.D0 - DABS(LONG)
        END IF
C       ----------------------------------------------------------------
        COLAT(K) = COLAT(K)*PI/180.D0
        ELONG(K) = ELONG(K)*PI/180.D0
            Z(K) = H*1.D-3
C       ----------------------------------------------------------------
       END IF
       K = K + 1
      END DO 
      KMAX = K - 1
C-----------------------------------------------------------------------
      JD = JDI
      PRINT *,'CALCULANDO Y ESCRIBIENDO ARCHIVO DE SALIDA.'
      DO WHILE (JD.LE.JDF)
C      -----------------------------------------------------------------
C      PROCESAMIENTO DE DATOS
C      -----------------------------------------------------------------
          D = JD - 2451545.D0
          T = D/36525.D0
        THP = 350.891982443297D0
       PROT = 360.D0/THP
         DT = FDT*PROT
C      -----------------------------------------------------------------
C      CALCULO DEL VECTOR DE POSICION DE MARTE RESPECTO AL SOL EN EL 
C      SISTEMA ICRF
C      -----------------------------------------------------------------
       CALL MARSORB(T,RME,I,OM)
C      -----------------------------------------------------------------
C      CALCULO DE LOS ANGULOS DE EULER QUE DAN LA ORIENTACION DE LA 
C      FIGURA DE MARTE EN EL SISTEMA ICRF
C      -----------------------------------------------------------------
        CALL MARSROT(D,T,ZETA,XI,THETA)
C
        XPNI =  DSIN(XI)*DSIN(ZETA)
        YPNI = -DSIN(XI)*DCOS(ZETA)
        ZPNI =  DCOS(XI)
C
        XPNE =  XPNI
        YPNE =  YPNI*DCOS(EPSE) + ZPNI*DSIN(EPSE)
        ZPNE = -YPNI*DSIN(EPSE) + ZPNI*DCOS(EPSE)
C
        XPNEI =  XPNE*DCOS(OM) + YPNE*DSIN(OM)
        YPNEI = -XPNE*DSIN(OM) + YPNE*DCOS(OM)
        ZPNEI =  ZPNE
C
        XPNO =  XPNEI
        YPNO =  YPNEI*DCOS(I) + ZPNEI*DSIN(I)
        ZPNO = -YPNEI*DSIN(I) + ZPNEI*DCOS(I)
C     
        LAM = DATAN2(XPNO,-YPNO)
C        EPS = DATAN2(PYTHAG(XPNO,YPNO),ZPNO)
C      -----------------------------------------------------------------
       IF (LROT1) THEN
        OBL = OBL0*PI/180.D0
C       ----------------------------------------------------------------
C       VERSOR MOMENTO ANGULAR ROTACIONAL EXPRESADO EN EL SISTEMA 
C       ECUATORIAL AREOCENTRICO.
C       ----------------------------------------------------------------
        HROTEA(1) =  DSIN(LAM)*DSIN(OBL)
        HROTEA(2) = -DCOS(LAM)*DSIN(OBL)
        HROTEA(3) =  DCOS(OBL)
C       ----------------------------------------------------------------
C       TRANSFORMACION A ICRF
C       ----------------------------------------------------------------
        HROTI(1) = HROTEA(1)
        HROTI(2) = HROTEA(2)*DCOS(I) - HROTEA(3)*DSIN(I)
        HROTI(3) = HROTEA(2)*DSIN(I) + HROTEA(3)*DCOS(I)
C       ----------------------------------------------------------------
        HROTO(1) = HROTI(1)*DCOS(OM) - HROTI(2)*DSIN(OM)
        HROTO(2) = HROTI(1)*DSIN(OM) + HROTI(2)*DCOS(OM)
        HROTO(3) = HROTI(3)
C       ----------------------------------------------------------------
        HROT(1) = HROTO(1)
        HROT(2) = HROTO(2)*DCOS(EPSE) - HROTO(3)*DSIN(EPSE)
        HROT(3) = HROTO(2)*DSIN(EPSE) + HROTO(3)*DCOS(EPSE)
C       ----------------------------------------------------------------
        ZETA = DATAN2(HROT(1),-HROT(2))
          XI = DATAN2(PYTHAG(HROT(1),HROT(2)),HROT(3))
C       ----------------------------------------------------------------
       ELSE IF (LROT2) THEN
        CALL MARSROT(D,T,ZETA,XI,THETA)
       ELSE
        PRINT *,'ERROR EN EL NUMERO IDENTIFICADOR DEL MODELO ROTACIONAL'
       END IF
C      -----------------------------------------------------------------
C      TRANFORMACION A SISTEMA ECUATORIAL AREOCENTRICO
C      -----------------------------------------------------------------
       RMME(1) =         RME(1)*DCOS(ZETA) + RME(2)*DSIN(ZETA)
       RMME(2) = (-1.D0)*RME(1)*DSIN(ZETA) + RME(2)*DCOS(ZETA)
       RMME(3) =         RME(3)
C
       RMMZ(1) =         RMME(1)
       RMMZ(2) =         RMME(2)*DCOS(XI) + RMME(3)*DSIN(XI)
       RMMZ(3) = (-1.D0)*RMME(2)*DSIN(XI) + RMME(3)*DCOS(XI)
C
        RMM(1) =         RMMZ(1)*DCOS(THETA) + RMMZ(2)*DSIN(THETA)
        RMM(2) = (-1.D0)*RMMZ(1)*DSIN(THETA) + RMMZ(2)*DCOS(THETA)
        RMM(3) =         RMMZ(3)
C-----------------------------------------------------------------------
C      ARCHIVOS DE SALIDA
C-----------------------------------------------------------------------
       DATENUM = JD-JD0-1.D0
       WRITE(CHDIA,'(I6)') INT(DATENUM)
       WRITE(CHHORA,'(I5)') INT((DATENUM - INT(DATENUM))*1D5)
C      -----------------------------------------------------------------
       IF (TRIM(PATHOUT).EQ.'.') THEN
        OPEN(UNIT=12,FILE='MARTESOL_'//TRIM(CHDIA)//'-'//TRIM(CHHORA)//
     &       '.OUT',STATUS='UNKNOWN')
       ELSE 
        OPEN(UNIT=12,FILE=TRIM(PATHOUT)//'MARTESOL_'//TRIM(CHDIA)//'-'
     &       //TRIM(CHHORA)//'.OUT',STATUS='UNKNOWN')
       END IF    
10     FORMAT (' ',F12.7,',',1X,F12.7,',',1X,F10.5,',',1X,F10.5,',',1X,
     &             F11.9)
C-----------------------------------------------------------------------
       DO K=1,KMAX
C       ----------------------------------------------------------------
C       DETERMINACIÓN DE LAS COMPONENTES DEL VECTOR R
C       (PUNTO SOBRE LA SUPERFICIE DE MARTE)
C       ----------------------------------------------------------------
        R(1) = (RADM+Z(K))*DSIN(COLAT(K))*DCOS(ELONG(K))/UAKM
        R(2) = (RADM+Z(K))*DSIN(COLAT(K))*DSIN(ELONG(K))/UAKM
        R(3) = (RADM+Z(K))*DCOS(COLAT(K))/UAKM
C       ----------------------------------------------------------------
C       DETERMINACION DE LAS COMPONENTES DEL VECTOR R_S
C       (POSICION DEL SOL EN EL PUNTO ELEGIDO DE LA SUPERFICIE DE MARTE)
C       ----------------------------------------------------------------
        RS(1) = -(RMM(1) + R(1))
        RS(2) = -(RMM(2) + R(2))
        RS(3) = -(RMM(3) + R(3))
          DMS = PYTHAG(PYTHAG(RS(1),RS(2)),RS(3))
C       ----------------------------------------------------------------
C       TRANSFORMACIÓN A SISTEMA HORIZONTAL
C       ----------------------------------------------------------------
        RSLG(1) =  RS(1)*DCOS(ELONG(K)) + RS(2)*DSIN(ELONG(K))
        RSLG(2) = -RS(1)*DSIN(ELONG(K)) + RS(2)*DCOS(ELONG(K))
        RSLG(3) =  RS(3)
C       
        RSLT(1) =  RSLG(1)*DCOS(COLAT(K)) - RSLG(3)*DSIN(COLAT(K))
        RSLT(2) =  RSLG(2)
        RSLT(3) =  RSLG(1)*DSIN(COLAT(K)) + RSLG(3)*DCOS(COLAT(K))
C
         RSH(1) =  RSLT(2)
         RSH(2) = -RSLT(1)
         RSH(3) =  RSLT(3)
C       ----------------------------------------------------------------
C       DETERMINACION DEL ACIMUT Y DISTANCIA CENITAL DEL SOL EN EL PUNTO
C       ELEGIDO DE LA SUPERFICIE DE MARTE
C       ----------------------------------------------------------------
C        IF ((RSH(1).GT.0.D0).AND.(RSH(2).GT.0.D0)) THEN
C            AC = DATAN(DABS(RSH(1)/RSH(2)))
C        ELSE IF ((RSH(1).GT.0.D0).AND.(RSH(2).LT.0.D0)) THEN
C            AC = PI - DATAN(DABS(RSH(1)/RSH(2)))
C        ELSE IF ((RSH(1).LT.0.D0).AND.(RSH(2).LT.0.D0)) THEN
C            AC = DATAN(DABS(RSH(1)/RSH(2))) + PI
C        ELSE IF ((RSH(1).LT.0.D0).AND.(RSH(2).GT.0.D0)) THEN
C            AC = 2.D0*PI - DATAN(DABS(RSH(1)/RSH(2)))
C        END IF
         AC = DATAN2(RSH(1),RSH(2))
         AC = AC*180.D0/PI
         IF (AC.LT.0.D0) AC = AC + 360.D0
        ALT = DATAN2(RSH(3),PYTHAG(RSH(1),RSH(2)))
        ALT = ALT*180.D0/PI
C       ----------------------------------------------------------------
C       ESCRITURA EN ARCHIVO DE SALIDA
C       ----------------------------------------------------------------
C                     1    2  3   4   5
        WRITE(12,10) LONG,LAT,AC,ALT,DMS
C       ----------------------------------------------------------------
       END DO
       CLOSE(UNIT=12)
       JD = JD + DT
      END DO
      PRINT *, 'FINALIZADO.'
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE CAL2JD(DD,MM,YY,HR,MIN,SEGS,JD)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER DD,MM,YY,HR,MIN,YD,MD,A,B,C,D
      REAL*8 SEGS,DR,TEST,GREG,JD
C-----------------------------------------------------------------------
      DR = DBLE(DD) + (DBLE(HR) + (DBLE(MIN) + SEGS/60.D0)/60.D0)/24.D0
      IF (MM.LT.3) THEN
            YD = YY - 1
            MD = MM + 12
      ELSE
            YD = YY
            MD = MM
      END IF
      A = INT(DBLE(YD)/100.D0)
      GREG = 1582.D0 + 10.D0/12.D0 + 15.0/365.25D0
      TEST = DBLE(YY) + DBLE(MD)/12.D0 + DR/365.25D0
      IF (TEST.GT.GREG) THEN
            B = 2 - A + INT(DBLE(A)/4.D0)
      ELSE
            B = 0
      END IF
      IF (YD.LT.0) THEN
            C = INT(365.25D0*DBLE(YD)-0.75D0)
      ELSE 
            C = INT(365.25D0*DBLE(YD))
      END IF
      D = INT(30.6001D0*DBLE(MD+1))
      JD = B+C+D+DR+1720994.5D0
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE MARSORB(T,RME,I,OM)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      REAL*8 PI,T,RM(3),EPSE,A0,DADT,E0,DEDT,I0,DIDT,L0,DLDT,VW0,DWDT,
     &       OM0,DOMDT,A,E,I,L,VW,OM,W,AM,PRIVUE,AE,RR,XO,YO,RME(3)
      PARAMETER ( PI = 4.D0*DATAN(1.D0) )
C-----------------------------------------------------------------------
C     OBLICUIDAD DE LA TIERRA
C-----------------------------------------------------------------------
      EPSE = 23.43928D0*PI/180.D0
C-----------------------------------------------------------------------
C     PARAMETROS ORBITALES DE MARTE
C-----------------------------------------------------------------------
         A0 =  1.52371034D0
       DADT =  0.00001847D0
         E0 =  0.09339410D0
       DEDT =  0.00007882D0
         I0 =  1.84969142D0
       DIDT = -0.00813131D0
         L0 = -4.55343205D0
       DLDT = 19140.30268499D0
        VW0 = -23.94362959D0
       DWDT = 0.44441088D0
        OM0 = 49.55953891D0
      DOMDT = -0.29257343D0
C     ------------------------------------------------------------------
       A =  A0 + DADT*T
       E =  E0 + DEDT*T
       I =  I0 + DIDT*T
       L =  L0 + DLDT*T
      VW = VW0 + DWDT*T
      OM = OM0 + DOMDT*T
       W = VW - OM
      AM = L - VW
C-----------------------------------------------------------------------
C     CONVERSION A RADIANES
C-----------------------------------------------------------------------
       I = PRIVUE(I)*PI/180.D0
       W = PRIVUE(W)*PI/180.D0
      OM = PRIVUE(OM)*PI/180.D0
      AM = PRIVUE(AM)*PI/180.D0
C-----------------------------------------------------------------------
C     CALCULO DE LA ANOMALIA VERDADERA
C-----------------------------------------------------------------------
      CALL SOLKEP(E,AM,AE)
C-----------------------------------------------------------------------
C     CALCULO DE LA DISTANCIA MARTE-SOL
C-----------------------------------------------------------------------
      RR = A*(1.D0-E*DCOS(AE))
C-----------------------------------------------------------------------
C     DETERMINACION DE LAS COMPONENTES DEL VECTOR R_M
C     (POSICION DEL SOL RESPECTO A MARTE EN SISTEMA ECLIPTICO)
C-----------------------------------------------------------------------
      XO = A*(DCOS(AE)-E)
      YO = A*DSQRT(1.D0-E*E)*DSIN(AE)
      RM(1) = (DCOS(OM)*DCOS(W)-DSIN(OM)*DSIN(W)*DCOS(I))*XO -
     &        (DSIN(W)*DCOS(OM)+DCOS(W)*DSIN(OM)*DCOS(I))*YO
      RM(2) = (DCOS(W)*DSIN(OM)+DSIN(W)*DCOS(OM)*DCOS(I))*XO +
     &        (-DSIN(W)*DSIN(OM)+DCOS(W)*DCOS(OM)*DCOS(I))*YO
      RM(3) = DSIN(W)*DSIN(I)*XO + DCOS(W)*DSIN(I)*YO
C-----------------------------------------------------------------------
C     TRANSFORMACION A ICRF
C-----------------------------------------------------------------------
      RME(1) = RM(1)
      RME(2) = RM(2)*DCOS(EPSE) - RM(3)*DSIN(EPSE)
      RME(3) = RM(2)*DSIN(EPSE) + RM(3)*DCOS(EPSE)
C-----------------------------------------------------------------------
      RETURN
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE MARSROT(D,T,ZETA,EPS,THETA)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      REAL*8 PI,D,T,ALPHA0,DELTA0,THETA,PRIVUE,ZETA,EPS
      PARAMETER ( PI = 4.D0*DATAN(1.D0) )
C-----------------------------------------------------------------------
C     PARAMETROS ROTACIONALES DE MARTE
C-----------------------------------------------------------------------
      ALPHA0 = 317.269202D0 - 0.10927547D0*T 
     &       + 0.000068D0*DSIN(198.991226D0 + 19139.4819985D0*T)
     &       + 0.000238D0*DSIN(226.292679D0 + 38280.8511281D0*T)
     &       + 0.000052D0*DSIN(249.663391D0 + 57420.7251593D0*T)
     &       + 0.000009D0*DSIN(266.183510D0 + 76560.6367950D0*T)
     &       + 0.419057D0*DSIN(79.398797D0 + 0.5042615D0*T)
      DELTA0 = 54.432516D0 - 0.05827105*T 
     &       + 0.000051D0*DCOS(122.433576D0 + 19139.9407476D0*T)
     &       + 0.000141D0*DCOS(43.058401D0 + 38280.8753272D0*T)
     &       + 0.000031D0*DCOS(57.663379D0 + 57420.7517205*T)
     &       + 0.000005D0*DCOS(79.476401D0 + 76560.6495004*T)
     &       + 1.591274D0*DCOS(166.325722D0 + 0.5042615*T)
       THETA = 176.049863D0 + 350.891982443297*D
     &       + 0.000145D0*DSIN(129.071773D0 + 19140.0328244D0*T)
     &       + 0.000157D0*DSIN(36.352167D0 + 38281.0473591D0*T)
     &       + 0.000040D0*DSIN(56.668646D0 + 57420.9295360D0*T)
     &       + 0.000001D0*DSIN(67.364003D0 + 76560.2552215D0*T)
     &       + 0.000001D0*DSIN(104.792680D0 + 95700.4387578D0*T)
     &       + 0.584542D0*DSIN(95.391654D0 + 0.5042615*T)
C     ------------------------------------------------------------------
      ALPHA0 = PRIVUE(ALPHA0)
      DELTA0 = PRIVUE(DELTA0)
C     ------------------------------------------------------------------
        ZETA = ALPHA0 + 90.D0
      IF (DELTA0.GE.0.D0) THEN
         EPS = 90.D0 - DELTA0
      ELSE
         EPS = 90.D0 + DABS(DELTA0)
      END IF
C     ------------------------------------------------------------------
       ZETA = PRIVUE(ZETA)*PI/180.D0
        EPS = PRIVUE(EPS)*PI/180.D0
      THETA = PRIVUE(THETA)*PI/180.D0
C-----------------------------------------------------------------------
      RETURN
C-----------------------------------------------------------------------
      END
C=======================================================================
      FUNCTION PYTHAG(A,B)
      REAL*8 A,B,PYTHAG
      REAL*8 ABSA,ABSB
      ABSA=ABS(A)
      ABSB=ABS(B)
      IF(ABSA.GT.ABSB)THEN
        PYTHAG=ABSA*SQRT(1.D0+(ABSB/ABSA)**2)
      ELSE
        IF(ABSB.EQ.0.D0)THEN
          PYTHAG=0.D0
        ELSE
          PYTHAG=ABSB*SQRT(1.D0+(ABSA/ABSB)**2)
        ENDIF
      ENDIF
      RETURN
      END
C=======================================================================
      FUNCTION PRIVUE(X)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      REAL*8 PI,X,PRIVUE
      PARAMETER (PI = 4.D0*DATAN(1.D0))
C-----------------------------------------------------------------------
      IF (X.LE.0.D0) X = X + 360.D0
      IF (X.GE.360.D0) THEN
       PRIVUE = X - DBLE(INT(X/360.D0))*360.D0
      ELSE
       PRIVUE = X
      END IF
C-----------------------------------------------------------------------
      RETURN
C-----------------------------------------------------------------------
      END
C=======================================================================
C SOLUCION DE LA ECUACION DE KEPLER ELIPTICA:
      SUBROUTINE SOLKEP(EX,M,E)
C
C SOLUCION ITERATIVA DE LA ECUACION DE KEPLER
C ENTRA:EX   EXCENTRICIDAD            (<1)
C       M    ANOMALIA MEDIA           (RADIANES)
C SALE: E    ANOMALIA EXCENTRICA      (RADIANES)
C
       IMPLICIT REAL*8 (A-H,O-Z)
       REAL*8 M,MK

       TOLE=1.D-10
       DPI=8.D0*DATAN(1.D0)
       M=DMOD(M,DPI)
       E=M
       NITER=0
 
 100     E0=E
      
       SE=DSIN(E0)
       CE=DCOS(E0)
       ES=EX*SE
       EC=1.D0-EX*CE
       MK=E0-ES
       U=(MK-M)/EC
       XPRI=E0-U
       XSEG=E0-U/(1.D0-U*ES)
       E=(XPRI+XSEG)/2.D0
       DEX=DABS(E-E0)
      
       NITER=NITER+1
       IF(NITER.GT.20)GOTO 200
       IF(DEX.GT.TOLE)GOTO 100                
       RETURN

C SI EL NUMERO DE ITERACIONES ES > 20 PRUEBA CON BISECCION:

C METODO DICOTOMICO:
200      CONTINUE        
       NDIC=0
       E0=-DPI
       DE0=DPI/10.D0

400     DE=DE0/(10.D0**NDIC)
       SE=DSIN(E0)
       CE=DCOS(E0)
       ES=EX*SE
       EM0=E0-ES-M
       NITER=0
      
300      E1=E0+DE
       NITER=NITER+1
      
       IF(NITER.GT.100)THEN
C      WRITE(50,*)'ERROR EN LA SOLUCION DE LA ECUACION DE KEPLER'
       RETURN
       ENDIF

       SE=DSIN(E1)
       CE=DCOS(E1)
       ES=EX*SE
       EM1=E1-ES-M
       IF(EM1*EM0.GT.0.D0)THEN
       E0=E1
       EM0=EM1
       GOTO 300
       ELSE
       NDIC=NDIC+1
       IF(NDIC.EQ.3)THEN
       E=E1
       RETURN
       ENDIF
       GOTO 400
      
      ENDIF

       RETURN
      END
C=======================================================================

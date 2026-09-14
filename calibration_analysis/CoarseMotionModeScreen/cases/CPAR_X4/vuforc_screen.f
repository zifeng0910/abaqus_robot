C=======================================================================
!DEC$ OBJCOMMENT LIB:'J:\abaqusfangzhen\socket_win32.lib'
C Reduced-Hydro VUAMP bridge.
C Magnetic SOCKET_* amplitudes are unchanged from production.  Six
C HYDRO_* amplitudes are calculated locally from RP V/VR sensors using a
C body-following, non-negative linear drag/rotational damping model.
C Consistent units: mm, tonne, s, N, N*mm.  Hydro power is non-positive.
C=======================================================================
      SUBROUTINE VUAMP(ampName,time,ampValueOld,dt,nprops,props,
     * nSvars,svars,lFlagsInfo,nSensor,sensorValues,sensorNames,
     * jSensorLookUpTable,AmpValueNew,lFlagsDefine,AmpDerivative,
     * AmpSecDerivative,AmpIncIntegral)
      INCLUDE 'VABA_PARAM.INC'
      PARAMETER (iStepTime=1,iTotalTime=2,nTime=2)
      PARAMETER (iInitialization=1,iRegularInc=2,ikStep=3,nFlagsInfo=3)
      PARAMETER (iComputeDeriv=1,iComputeSecDeriv=2,iComputeInteg=3,
     * iStopAnalysis=4,iConcludeStep=5,nFlagsDefine=5)
      DIMENSION time(nTime),lFlagsInfo(nFlagsInfo),
     * lFlagsDefine(nFlagsDefine),sensorValues(nSensor),props(nprops),
     * svars(nSvars),jSensorLookUpTable(*)
      CHARACTER*80 ampName,sensorNames(nSensor)
      DOUBLE PRECISION loads(6),hloads(6),lastTime,tSend
      DOUBLE PRECISION u1,u2,u3,r1,r2,r3,v1,v2,v3,w1,w2,w3
      DOUBLE PRECISION a0(3),a(3),v(3),w(3),vp,vs(3),wws(3),www(3)
      DOUBLE PRECISION cpar,cperp,kspin,kwob,ang,ss,cc,kx,ky,kz
      INTEGER connected,failCount,rc,idx,isMaster,initialized
      INTEGER logInit
      COMMON /SOCKBRIDGE/ loads,hloads,lastTime,connected,failCount
      SAVE /SOCKBRIDGE/
      SAVE initialized,logInit
      DATA initialized /0/, logInit /0/
      DOUBLE PRECISION VGETSENSORVALUE
      INTEGER SOCKET_OPEN,SOCKET_POSE
      EXTERNAL SOCKET_OPEN,SOCKET_POSE,VGETSENSORVALUE
C     Provisional coefficients calibrated to CSF scale; all non-negative.
C     C: N*s/mm. K: N*mm*s.  They are deliberately dissipative.
      DATA cpar /1.6e-08/, cperp /1.2D-8/
      DATA kspin /1.0D-9/, kwob /3e-09/
C     Directed HEAD->TAIL axis from the actual robot exterior geometry.
      DATA a0 /0.9647382600216D0,-0.1188742372140D0,
     *          0.2348382536499D0/
      DO 5 idx=1,nFlagsDefine
         lFlagsDefine(idx)=0
    5 CONTINUE
      AmpDerivative=0.0D0
      AmpSecDerivative=0.0D0
      AmpIncIntegral=0.0D0
      IF (initialized.EQ.0) THEN
         lastTime=-1.0D99
         connected=0
         failCount=0
         DO 10 idx=1,6
            loads(idx)=0.0D0
            hloads(idx)=0.0D0
   10    CONTINUE
         initialized=1
      ENDIF
      isMaster=1
      IF (lFlagsInfo(iInitialization).EQ.1) isMaster=0
      tSend=DBLE(time(iTotalTime))
      IF (isMaster.EQ.1 .AND.
     *    ABS(tSend-lastTime).GT.1.0D-15) THEN
C        Always fetch by the named sensor table.  This avoids relying on
C        implementation-dependent ordering of sensorValues(1:12).
         u1=VGETSENSORVALUE('RP_U1',jSensorLookUpTable,sensorValues)
         u2=VGETSENSORVALUE('RP_U2',jSensorLookUpTable,sensorValues)
         u3=VGETSENSORVALUE('RP_U3',jSensorLookUpTable,sensorValues)
         r1=VGETSENSORVALUE('RP_UR1',jSensorLookUpTable,sensorValues)
         r2=VGETSENSORVALUE('RP_UR2',jSensorLookUpTable,sensorValues)
         r3=VGETSENSORVALUE('RP_UR3',jSensorLookUpTable,sensorValues)
         v1=VGETSENSORVALUE('RP_V1',jSensorLookUpTable,sensorValues)
         v2=VGETSENSORVALUE('RP_V2',jSensorLookUpTable,sensorValues)
         v3=VGETSENSORVALUE('RP_V3',jSensorLookUpTable,sensorValues)
         w1=VGETSENSORVALUE('RP_VR1',jSensorLookUpTable,sensorValues)
         w2=VGETSENSORVALUE('RP_VR2',jSensorLookUpTable,sensorValues)
         w3=VGETSENSORVALUE('RP_VR3',jSensorLookUpTable,sensorValues)
C        Magnetic query remains one production Socket request/increment.
         IF (connected.EQ.0) THEN
            rc=SOCKET_OPEN()
            IF (rc.EQ.0) connected=1
         ENDIF
         IF (connected.EQ.1) THEN
            rc=SOCKET_POSE(tSend,-7.468174204284D0+u1,
     *       -3.676918015967D0+u2,-9.550745259298D0+u3,r1,r2,r3,loads)
            IF (rc.NE.0) THEN
               connected=0; failCount=failCount+1
               DO 20 idx=1,6
                  loads(idx)=0.0D0
   20          CONTINUE
            ELSE
               failCount=0
            ENDIF
         ELSE
            failCount=failCount+1
         ENDIF
C        Rodrigues rotation of initial body axis by Abaqus UR vector.
         ang=DSQRT(r1*r1+r2*r2+r3*r3)
         IF (ang.LT.1.0D-14) THEN
            a(1)=a0(1); a(2)=a0(2); a(3)=a0(3)
         ELSE
            kx=r1/ang; ky=r2/ang; kz=r3/ang
            cc=DCOS(ang); ss=DSIN(ang)
            a(1)=(cc+kx*kx*(1.0D0-cc))*a0(1)+
     *           (kx*ky*(1.0D0-cc)-kz*ss)*a0(2)+
     *           (kx*kz*(1.0D0-cc)+ky*ss)*a0(3)
            a(2)=(ky*kx*(1.0D0-cc)+kz*ss)*a0(1)+
     *           (cc+ky*ky*(1.0D0-cc))*a0(2)+
     *           (ky*kz*(1.0D0-cc)-kx*ss)*a0(3)
            a(3)=(kz*kx*(1.0D0-cc)-ky*ss)*a0(1)+
     *           (kz*ky*(1.0D0-cc)+kx*ss)*a0(2)+
     *           (cc+kz*kz*(1.0D0-cc))*a0(3)
         ENDIF
         v(1)=v1; v(2)=v2; v(3)=v3
         w(1)=w1; w(2)=w2; w(3)=w3
         vp=v(1)*a(1)+v(2)*a(2)+v(3)*a(3)
         vs(1)=v(1)-vp*a(1); vs(2)=v(2)-vp*a(2)
         vs(3)=v(3)-vp*a(3)
         hloads(1)=-cpar*vp*a(1)-cperp*vs(1)
         hloads(2)=-cpar*vp*a(2)-cperp*vs(2)
         hloads(3)=-cpar*vp*a(3)-cperp*vs(3)
         wws(1)=(w(1)*a(1)+w(2)*a(2)+w(3)*a(3))*a(1)
         wws(2)=(w(1)*a(1)+w(2)*a(2)+w(3)*a(3))*a(2)
         wws(3)=(w(1)*a(1)+w(2)*a(2)+w(3)*a(3))*a(3)
         www(1)=w(1)-wws(1); www(2)=w(2)-wws(2)
         www(3)=w(3)-wws(3)
         hloads(4)=-kspin*wws(1)-kwob*www(1)
         hloads(5)=-kspin*wws(2)-kwob*www(2)
         hloads(6)=-kspin*wws(3)-kwob*www(3)
         IF (logInit.EQ.0) THEN
            OPEN(77,FILE='hydro_increment.csv',
     *           STATUS='REPLACE')
            WRITE(77,*) 'time_s,v1_mm_s,v2_mm_s,v3_mm_s,vr1,vr2,vr3,',
     *       'Fh1_N,Fh2_N,Fh3_N,Th1_Nmm,Th2_Nmm,Th3_Nmm'
            logInit=1
         ENDIF
         WRITE(77,100) tSend,v1,v2,v3,w1,w2,w3,
     *       hloads(1),hloads(2),hloads(3),hloads(4),hloads(5),hloads(6)
  100    FORMAT(ES16.8,12(',',ES16.8))
         lastTime=tSend
      ENDIF
      IF (isMaster.EQ.1 .AND. failCount.GE.100)
     *    lFlagsDefine(iStopAnalysis)=1
      AmpValueNew=0.0D0
      IF (ampName(1:9).EQ.'SOCKET_FX') AmpValueNew=loads(1)
      IF (ampName(1:9).EQ.'SOCKET_FY') AmpValueNew=loads(2)
      IF (ampName(1:9).EQ.'SOCKET_FZ') AmpValueNew=loads(3)
      IF (ampName(1:9).EQ.'SOCKET_MX') AmpValueNew=loads(4)
      IF (ampName(1:9).EQ.'SOCKET_MY') AmpValueNew=loads(5)
      IF (ampName(1:9).EQ.'SOCKET_MZ') AmpValueNew=loads(6)
      IF (ampName(1:8).EQ.'HYDRO_FX') AmpValueNew=hloads(1)
      IF (ampName(1:8).EQ.'HYDRO_FY') AmpValueNew=hloads(2)
      IF (ampName(1:8).EQ.'HYDRO_FZ') AmpValueNew=hloads(3)
      IF (ampName(1:8).EQ.'HYDRO_MX') AmpValueNew=hloads(4)
      IF (ampName(1:8).EQ.'HYDRO_MY') AmpValueNew=hloads(5)
      IF (ampName(1:8).EQ.'HYDRO_MZ') AmpValueNew=hloads(6)
      RETURN
      END
C=======================================================================
      SUBROUTINE VEXTERNALDB(lOp,i_Array,niArray,r_Array,nrArray)
      INCLUDE 'VABA_PARAM.INC'
      DIMENSION i_Array(niArray),r_Array(nrArray)
      EXTERNAL SOCKET_CLOSE
      IF (lOp.EQ.6) THEN
         CALL SOCKET_CLOSE()
         CLOSE(77)
      ENDIF
      RETURN
      END

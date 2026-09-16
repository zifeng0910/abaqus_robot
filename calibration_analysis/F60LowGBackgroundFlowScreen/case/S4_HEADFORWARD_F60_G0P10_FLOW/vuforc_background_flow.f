C=======================================================================
!DEC$ OBJCOMMENT LIB:'J:\abaqusfangzhen\socket_win32.lib'
C Reduced-Hydro VUAMP bridge with prescribed background flow for F60 screen.
C Magnetic SOCKET_* amplitudes are unchanged from production. Six
C HYDRO_* amplitudes use V_rel=V_robot-U_flow*c_hat for translation;
C rotational damping remains body-following and unchanged.
C Consistent units: mm, tonne, s, N, N*mm. This is diagnostic screening only.
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
      DOUBLE PRECISION a0(3),c0(3),flow0(3),a(3),v(3),vrel(3),w(3)
      DOUBLE PRECISION vp,vs(3),wws(3),www(3)
      DOUBLE PRECISION vrobot_s,vrel_s,fmag_s,fhydro_s
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
      DATA cpar /4e-09/, cperp /1.2D-8/
      DATA kspin /1.0D-9/, kwob /3e-09/
C     Directed HEAD->TAIL axis from the actual robot exterior geometry.
      DATA a0 /-0.9762799602464D0,0.0618320247447D0,
     *          -0.2074951564187D0/
C     Prescribed diagnostic background flow: +10 mm/s along canonical +s.
      DATA c0 /0.9762799602464D0,-0.0618320247447D0,
     *          0.2074951564187D0/
      DATA flow0 /9.762799602464D0,-0.618320247447D0,
     *            2.074951564187D0/
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
         vrel(1)=v(1)-flow0(1); vrel(2)=v(2)-flow0(2)
         vrel(3)=v(3)-flow0(3)
         vp=vrel(1)*a(1)+vrel(2)*a(2)+vrel(3)*a(3)
         vs(1)=vrel(1)-vp*a(1); vs(2)=vrel(2)-vp*a(2)
         vs(3)=vrel(3)-vp*a(3)
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
            OPEN(77,FILE='J:\\abaqusfangzhen\\abaqus_robot\\calibration_analysis\\F60LowG'//
     *           'BackgroundFlowScreen\\case\\S4_HEADFORWARD_F60_G0P10_FLOW\\'//
     *           'hydro_increment.csv',STATUS='REPLACE')
            WRITE(77,*) 'time_s,v1_mm_s,v2_mm_s,v3_mm_s,',
     *       'flow1_mm_s,flow2_mm_s,flow3_mm_s,',
     *       'vrel1_mm_s,vrel2_mm_s,vrel3_mm_s,vr1,vr2,vr3,',
     *       'Fmag1_N,Fmag2_N,Fmag3_N,Fh1_N,Fh2_N,Fh3_N,',
     *       'Th1_Nmm,Th2_Nmm,Th3_Nmm,vrobot_s_mm_s,',
     *       'vrel_s_mm_s,U_flow_s_mm_s,Fmag_s_N,Fhydro_s_N'
            logInit=1
         ENDIF
         vrobot_s=v1*c0(1)+v2*c0(2)+v3*c0(3)
         vrel_s=vrel(1)*c0(1)+vrel(2)*c0(2)+vrel(3)*c0(3)
         fmag_s=loads(1)*c0(1)+loads(2)*c0(2)+loads(3)*c0(3)
         fhydro_s=hloads(1)*c0(1)+hloads(2)*c0(2)+hloads(3)*c0(3)
         WRITE(77,100) tSend,v1,v2,v3,flow0(1),flow0(2),flow0(3),
     *       vrel(1),vrel(2),vrel(3),w1,w2,w3,
     *       loads(1),loads(2),loads(3),hloads(1),hloads(2),hloads(3),
     *       hloads(4),hloads(5),hloads(6),vrobot_s,vrel_s,10.0D0,
     *       fmag_s,fhydro_s
  100    FORMAT(ES16.8,27(',',ES16.8))
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

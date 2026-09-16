C=======================================================================
!DEC$ OBJCOMMENT LIB:'J:\abaqusfangzhen\socket_win32.lib'
C Pure Magpylib socket VUAMP bridge for the F60 true-CEL case.
C No ReducedHydro coefficients, amplitudes, or external fluid loads exist.
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
      DOUBLE PRECISION loads(6),lastTime,tSend
      DOUBLE PRECISION u1,u2,u3,r1,r2,r3
      INTEGER connected,failCount,rc,idx,isMaster,initialized
      COMMON /SOCKBRIDGE/ loads,lastTime,connected,failCount
      SAVE /SOCKBRIDGE/
      SAVE initialized
      DATA initialized /0/
      DOUBLE PRECISION VGETSENSORVALUE
      INTEGER SOCKET_OPEN,SOCKET_POSE
      EXTERNAL SOCKET_OPEN,SOCKET_POSE,VGETSENSORVALUE
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
   10    CONTINUE
         initialized=1
      ENDIF
      isMaster=1
      IF (lFlagsInfo(iInitialization).EQ.1) isMaster=0
      tSend=DBLE(time(iTotalTime))
      IF (isMaster.EQ.1 .AND.
     *    ABS(tSend-lastTime).GT.1.0D-15) THEN
         u1=VGETSENSORVALUE('RP_U1',jSensorLookUpTable,sensorValues)
         u2=VGETSENSORVALUE('RP_U2',jSensorLookUpTable,sensorValues)
         u3=VGETSENSORVALUE('RP_U3',jSensorLookUpTable,sensorValues)
         r1=VGETSENSORVALUE('RP_UR1',jSensorLookUpTable,sensorValues)
         r2=VGETSENSORVALUE('RP_UR2',jSensorLookUpTable,sensorValues)
         r3=VGETSENSORVALUE('RP_UR3',jSensorLookUpTable,sensorValues)
         IF (connected.EQ.0) THEN
            rc=SOCKET_OPEN()
            IF (rc.EQ.0) connected=1
         ENDIF
         IF (connected.EQ.1) THEN
            rc=SOCKET_POSE(tSend,-7.468174204284D0+u1,
     *       -3.676918015967D0+u2,-9.550745259298D0+u3,r1,r2,r3,loads)
            IF (rc.NE.0) THEN
               connected=0
               failCount=failCount+1
               DO 20 idx=1,6
                  loads(idx)=0.0D0
   20          CONTINUE
            ELSE
               failCount=0
            ENDIF
         ELSE
            failCount=failCount+1
         ENDIF
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
      RETURN
      END
C=======================================================================
      SUBROUTINE VEXTERNALDB(lOp,i_Array,niArray,r_Array,nrArray)
      INCLUDE 'VABA_PARAM.INC'
      DIMENSION i_Array(niArray),r_Array(nrArray)
      EXTERNAL SOCKET_CLOSE
      IF (lOp.EQ.6) CALL SOCKET_CLOSE()
      RETURN
      END

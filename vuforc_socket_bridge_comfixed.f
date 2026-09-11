C=======================================================================
!DEC$ OBJCOMMENT LIB:'J:\abaqusfangzhen\socket_win32.lib'
C Historical file name retained for the requested workflow.  Abaqus/Explicit
C does NOT provide a VUFORC subroutine. This file correctly implements VUAMP
C (user amplitude for the six RP *CLOADs) and VEXTERNALDB (safe close).
C It obtains RP state through the six named *SENSOR output requests in the
C Socket-ready input deck. Sensor values are end-of-PREVIOUS-increment values,
C so the coupling is explicit staggered (one increment lag), not iterative.
C=======================================================================
      SUBROUTINE VUAMP(ampName,time,ampValueOld,dt,nprops,props,
     * nSvars,svars,lFlagsInfo,nSensor,sensorValues,sensorNames,
     * jSensorLookUpTable,AmpValueNew,lFlagsDefine,AmpDerivative,
     * AmpSecDerivative,AmpIncIntegral)
      INCLUDE 'VABA_PARAM.INC'
      PARAMETER (iStepTime=1,iTotalTime=2,nTime=2)
      PARAMETER (iInitialization=1,iRegularInc=2,ikStep=3,
     * nFlagsInfo=3)
      PARAMETER (iComputeDeriv=1,iComputeSecDeriv=2,
     * iComputeInteg=3,iStopAnalysis=4,iConcludeStep=5,
     * nFlagsDefine=5)
      DIMENSION time(nTime),lFlagsInfo(nFlagsInfo),
     * lFlagsDefine(nFlagsDefine),sensorValues(nSensor),
     * props(nprops),svars(nSvars),jSensorLookUpTable(*)
      CHARACTER*80 ampName,sensorNames(nSensor)
      DOUBLE PRECISION loads(6),lastTime,tSend,u1,u2,u3,r1,r2,r3
      INTEGER connected,failCount,rc,idx,initialized,isMaster
      COMMON /SOCKBRIDGE/ loads,lastTime,connected,failCount
      SAVE /SOCKBRIDGE/
      SAVE initialized
      DATA initialized /0/
      DOUBLE PRECISION VGETSENSORVALUE
      INTEGER SOCKET_OPEN,SOCKET_POSE
      EXTERNAL SOCKET_OPEN,SOCKET_POSE,VGETSENSORVALUE
C     Request no derivative/integral work from Abaqus for force amplitudes.
      DO 5 idx=1,nFlagsDefine
         lFlagsDefine(idx)=0
    5 CONTINUE
      AmpDerivative=0.0D0
      AmpSecDerivative=0.0D0
      AmpIncIntegral=0.0D0
C     Initialize the process-wide cache exactly once.  Abaqus calls VUAMP
C     separately for all six named amplitudes; using iInitialization here
C     reset the shared socket/cache up to six times.
      IF (initialized.EQ.0) THEN
         lastTime=-1.0D99
         connected=0
         failCount=0
         DO 10 idx=1,6
            loads(idx)=0.0D0
   10    CONTINUE
         initialized=1
      ENDIF
C     Only the first defined amplitude owns socket I/O.  The other five
C     amplitudes read the same cached six-component load vector.  This avoids
C     re-entrant request/response traffic during field-output callbacks.
      isMaster=1
C     Do not query the external model during VUAMP initialization: sensor
C     values are not valid until the first regular increment.
      IF (lFlagsInfo(iInitialization).EQ.1) isMaster=0
C     Abaqus may invoke the six user amplitudes in an implementation-
C     dependent order and can alter the callback name padding. Any
C     SOCKET_* amplitude may therefore own the de-duplicated request; the
C     lastTime guard below still guarantees one pose query per increment.
C     Query once per physical increment and cache six values for its six
C     VUAMP invocations.  VUAMP supplies state at increment start.
      tSend=DBLE(time(iTotalTime))
      IF (isMaster.EQ.1 .AND.
     *    ABS(tSend-lastTime).GT.1.0D-15) THEN
         IF (ABS(tSend).GT.10.0D0) THEN
C           Guard against a corrupted VUAMP time value.  Never send it to
C           the external field model; fail closed with zero loads.
            failCount=failCount+1
            DO 15 idx=1,6
               loads(idx)=0.0D0
   15       CONTINUE
            lastTime=tSend
            GOTO 30
         ENDIF
         IF (connected.EQ.0) THEN
            rc=SOCKET_OPEN()
            IF (rc.EQ.0) connected=1
         ENDIF
         IF (connected.EQ.1) THEN
C     The input deck defines the six RP sensor requests in a
C     fixed order (U1,U2,U3,UR1,UR2,UR3).  Read the values
C     directly; this avoids a name-table ambiguity seen in
C     the Windows VUAMP utility while retaining the exact
C     beginning-of-increment sensor semantics.
            IF (nSensor.GE.6) THEN
               u1=DBLE(sensorValues(1))
               u2=DBLE(sensorValues(2))
               u3=DBLE(sensorValues(3))
               r1=DBLE(sensorValues(4))
               r2=DBLE(sensorValues(5))
               r3=DBLE(sensorValues(6))
            ELSE
               u1=VGETSENSORVALUE('RP_U1',jSensorLookUpTable,
     *                            sensorValues)
               u2=VGETSENSORVALUE('RP_U2',jSensorLookUpTable,
     *                            sensorValues)
               u3=VGETSENSORVALUE('RP_U3',jSensorLookUpTable,
     *                            sensorValues)
               r1=VGETSENSORVALUE('RP_UR1',jSensorLookUpTable,
     *                            sensorValues)
               r2=VGETSENSORVALUE('RP_UR2',jSensorLookUpTable,
     *                            sensorValues)
               r3=VGETSENSORVALUE('RP_UR3',jSensorLookUpTable,
     *                            sensorValues)
            ENDIF
C           Send the RP in the active deck's ORIGINAL ABAQUS GLOBAL frame.
C           The Magpylib server is solely responsible for reading the one
C           shared R,t JSON, mapping pose to the flat-DXF frame, and mapping
C           returned force/torque back with R-transpose.  Do not duplicate
C           R,t constants in this Fortran bridge.
C           All lengths are mm; returned force and torque are N and N*mm in
C           the Abaqus global frame.
C           COM-fixed calibration deck RP_ROBOT at t=0.  This is the
C           actual Abaqus assembly coordinate; the server applies the shared
C           R,t transform to the flat-DXF/Magpylib frame.
            rc=SOCKET_POSE(tSend,-7.468174204284D0+u1,
     *           -3.676918015967D0+u2,-9.550745259298D0+u3,
     *           r1,r2,r3,loads)
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
   30 CONTINUE
C     A finite fail threshold prevents indefinite deadlock / stale forcing.
C     Tolerate transient socket delays; fail-safe loads are zeroed above.
C     Stop only after a prolonged outage rather than five 1-second stalls.
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
C     lOp=6 is Abaqus/Explicit end of analysis. Close TCP gracefully.
      IF (lOp.EQ.6) CALL SOCKET_CLOSE()
      RETURN
      END

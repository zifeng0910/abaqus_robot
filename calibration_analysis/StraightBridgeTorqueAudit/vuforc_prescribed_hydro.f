C=======================================================================
C Reduced-Hydro only VUAMP for the prescribed-rotation replay.
C=======================================================================
      SUBROUTINE VUAMP(ampName,time,ampValueOld,dt,nprops,props,
     * nSvars,svars,lFlagsInfo,nSensor,sensorValues,sensorNames,
     * jSensorLookUpTable,AmpValueNew,lFlagsDefine,AmpDerivative,
     * AmpSecDerivative,AmpIncIntegral)
      INCLUDE 'VABA_PARAM.INC'
      PARAMETER (iInitialization=1,iRegularInc=2,ikStep=3,nFlagsInfo=3)
      PARAMETER (iComputeDeriv=1,iComputeSecDeriv=2,iComputeInteg=3,
     * iStopAnalysis=4,iConcludeStep=5,nFlagsDefine=5)
      DIMENSION time(2),lFlagsInfo(nFlagsInfo),lFlagsDefine(nFlagsDefine),
     * props(nprops),svars(nSvars),jSensorLookUpTable(*)
      CHARACTER*80 ampName,sensorNames(nSensor)
      DOUBLE PRECISION h(6),u1,u2,u3,r1,r2,r3,v1,v2,v3,w1,w2,w3
      DOUBLE PRECISION a0(3),a(3),v(3),w(3),vp,vs(3),ws(3),ww(3)
      DOUBLE PRECISION cpar,cperp,kspin,kwob,ang,ss,cc,kx,ky,kz
      DOUBLE PRECISION VGETSENSORVALUE
      INTEGER i
      DATA cpar /4e-09/, cperp /1.2D-8/
      DATA kspin /1.0D-9/, kwob /3.0D-9/
      DATA a0 /0.9647382600216D0,-0.1188742372140D0,
     *          0.2348382536499D0/
      DO 5 i=1,nFlagsDefine
         lFlagsDefine(i)=0
    5 CONTINUE
      AmpDerivative=0.D0
      AmpSecDerivative=0.D0
      AmpIncIntegral=0.D0
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
      ang=DSQRT(r1*r1+r2*r2+r3*r3)
      IF (ang.LT.1.D-14) THEN
         a(1)=a0(1); a(2)=a0(2); a(3)=a0(3)
      ELSE
         kx=r1/ang; ky=r2/ang; kz=r3/ang
         cc=DCOS(ang); ss=DSIN(ang)
         a(1)=(cc+kx*kx*(1.D0-cc))*a0(1)+
     *    (kx*ky*(1.D0-cc)-kz*ss)*a0(2)+
     *    (kx*kz*(1.D0-cc)+ky*ss)*a0(3)
         a(2)=(ky*kx*(1.D0-cc)+kz*ss)*a0(1)+
     *    (cc+ky*ky*(1.D0-cc))*a0(2)+
     *    (ky*kz*(1.D0-cc)-kx*ss)*a0(3)
         a(3)=(kz*kx*(1.D0-cc)-ky*ss)*a0(1)+
     *    (kz*ky*(1.D0-cc)+kx*ss)*a0(2)+
     *    (cc+kz*kz*(1.D0-cc))*a0(3)
      ENDIF
      v(1)=v1; v(2)=v2; v(3)=v3
      w(1)=w1; w(2)=w2; w(3)=w3
      vp=v(1)*a(1)+v(2)*a(2)+v(3)*a(3)
      DO 10 i=1,3
         vs(i)=v(i)-vp*a(i)
         ws(i)=(w(1)*a(1)+w(2)*a(2)+w(3)*a(3))*a(i)
         ww(i)=w(i)-ws(i)
         h(i)=-cpar*vp*a(i)-cperp*vs(i)
         h(i+3)=-kspin*ws(i)-kwob*ww(i)
   10 CONTINUE
      AmpValueNew=0.D0
      IF (ampName(1:8).EQ.'HYDRO_FX') AmpValueNew=h(1)
      IF (ampName(1:8).EQ.'HYDRO_FY') AmpValueNew=h(2)
      IF (ampName(1:8).EQ.'HYDRO_FZ') AmpValueNew=h(3)
      IF (ampName(1:8).EQ.'HYDRO_MX') AmpValueNew=h(4)
      IF (ampName(1:8).EQ.'HYDRO_MY') AmpValueNew=h(5)
      IF (ampName(1:8).EQ.'HYDRO_MZ') AmpValueNew=h(6)
      RETURN
      END
C=======================================================================

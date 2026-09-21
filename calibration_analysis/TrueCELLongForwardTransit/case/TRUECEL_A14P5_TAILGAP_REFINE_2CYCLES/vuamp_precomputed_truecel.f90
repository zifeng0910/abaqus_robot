module truecel_precomputed_data
  implicit none
  integer, parameter :: max_s=2001, max_phase=721
  integer :: ns=0, nphase=0, lookup_count=0
  real(8) :: s_axis(max_s), phase_axis(max_phase)
  real(8) :: b_unit(max_phase,3), grad_unit(max_s,9)
  logical :: table_loaded=.false., log_open=.false.
contains
  subroutine load_magnetic_table()
    integer :: i,j,ios
    if (table_loaded) return
    open(unit=91,file=trim('J:\abaqusfangzhen\abaqus_robot\calibration_analysis\' // &
      'TrueCELLongForwardTransit\case\TRUECEL_A14P5_DIAG_2CYCLES\' // &
      'magnetic_field_gradient_table_A14P5.dat'), &
      status='old',action='read',iostat=ios)
    if (ios /= 0) then
      call xplb_abqerr(-3,'PRECOMPUTED_TABLE_FATAL: missing table',0,0.0d0,' ')
      return
    endif
    read(91,*,iostat=ios) ns,nphase
    if (ios /= 0 .or. ns<2 .or. ns>max_s .or. nphase<2 .or. nphase>max_phase) then
      call xplb_abqerr(-3,'PRECOMPUTED_TABLE_FATAL: dimensions',0,0.0d0,' ')
      return
    endif
    do i=1,nphase
      read(91,*,iostat=ios) phase_axis(i),(b_unit(i,j),j=1,3)
      if (ios /= 0) exit
    enddo
    do i=1,ns
      read(91,*,iostat=ios) s_axis(i),(grad_unit(i,j),j=1,9)
      if (ios /= 0) exit
    enddo
    close(91)
    if (ios /= 0) then
      call xplb_abqerr(-3,'PRECOMPUTED_TABLE_FATAL: truncated table',0,0.0d0,' ')
      return
    endif
    table_loaded=.true.
    write(6,*) 'PRECOMPUTED_TABLE_LOADED ns=',ns,' nphase=',nphase
  end subroutine load_magnetic_table

  subroutine interpolate_table(s_eff,phase_deg,b,grad,outside)
    real(8),intent(in) :: s_eff,phase_deg
    real(8),intent(out) :: b(3),grad(3,3)
    logical,intent(out) :: outside
    integer :: is,ip,i,j,k
    real(8) :: ws,wp,phase
    lookup_count=lookup_count+1
    outside=s_eff<s_axis(1) .or. s_eff>s_axis(ns)
    if (outside) then
      b=0.0d0; grad=0.0d0; return
    endif
    phase=modulo(phase_deg,360.0d0)
    is=min(ns-1,max(1,int((s_eff-s_axis(1))/(s_axis(2)-s_axis(1)))+1))
    ip=min(nphase-1,max(1,int((phase-phase_axis(1))/(phase_axis(2)-phase_axis(1)))+1))
    ws=(s_eff-s_axis(is))/(s_axis(is+1)-s_axis(is))
    wp=(phase-phase_axis(ip))/(phase_axis(ip+1)-phase_axis(ip))
    do i=1,3
      b(i)=(1.0d0-wp)*b_unit(ip,i)+wp*b_unit(ip+1,i)
    enddo
    k=0
    do i=1,3
      do j=1,3
        k=k+1
        grad(i,j)=(1.0d0-ws)*grad_unit(is,k)+ws*grad_unit(is+1,k)
      enddo
    enddo
  end subroutine interpolate_table

  subroutine open_increment_log()
    if (log_open) return
    open(unit=77,file='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\TrueCELLongForwardTransit\case\TRUECEL_A14P5_TAILGAP_REFINE_2CYCLES\magnetic_increment.csv',status='replace')
    write(77,'(A)') 'time_s,s_eff_mm,phase_deg,Fmag1_N,Fmag2_N,Fmag3_N,Mmag1_Nmm,Mmag2_Nmm,Mmag3_Nmm'
    log_open=.true.
  end subroutine open_increment_log
end module truecel_precomputed_data

subroutine vuamp(ampName,time,ampValueOld,dt,nprops,props,nSvars,svars, &
  lFlagsInfo,nSensor,sensorValues,sensorNames,jSensorLookUpTable, &
  AmpValueNew,lFlagsDefine,AmpDerivative,AmpSecDerivative,AmpIncIntegral)
  use truecel_precomputed_data
  include 'VABA_PARAM.INC'
  parameter (iStepTime=1,iTotalTime=2,nTime=2)
  parameter (iInitialization=1,iRegularInc=2,ikStep=3,nFlagsInfo=3)
  parameter (iComputeDeriv=1,iComputeSecDeriv=2,iComputeInteg=3,iStopAnalysis=4,iConcludeStep=5,nFlagsDefine=5)
  dimension time(nTime),lFlagsInfo(nFlagsInfo),lFlagsDefine(nFlagsDefine)
  dimension sensorValues(nSensor),props(nprops),svars(nSvars),jSensorLookUpTable(*)
  character*80 ampName,sensorNames(nSensor)
  real(8) :: loads(6),last_time,t,u(3),ur(3)
  real(8) :: c(3),a0(3),a(3),m0(3),moment(3)
  real(8) :: b(3),grad(3,3),ang,kx,ky,kz,cc,ss,q,scale,s_eff,phase_deg,dotu
  integer :: initialized,idx,i,j,is_master
  logical :: outside
  double precision VGETSENSORVALUE
  save loads,last_time,initialized
  data initialized /0/
  data c /0.9762799602464296d0,-0.0618320247446977d0,0.2074951564186524d0/
  data a0 /-0.9762799602464296d0,0.0618320247446977d0,-0.2074951564186524d0/
  data m0 /-1.061824299255672d-3,6.72499168508204d-5,-2.25676452654195d-4/
  do idx=1,nFlagsDefine
    lFlagsDefine(idx)=0
  enddo
  AmpDerivative=0.0d0; AmpSecDerivative=0.0d0; AmpIncIntegral=0.0d0
  if (.not.table_loaded) call load_magnetic_table()
  if (initialized==0) then
    loads=0.0d0; last_time=-1.0d99; initialized=1
  endif
  is_master=1
  if (lFlagsInfo(iInitialization)==1) is_master=0
  t=dble(time(iTotalTime))
  if (is_master==1 .and. abs(t-last_time)>1.0d-15) then
    u(1)=VGETSENSORVALUE('RP_U1',jSensorLookUpTable,sensorValues)
    u(2)=VGETSENSORVALUE('RP_U2',jSensorLookUpTable,sensorValues)
    u(3)=VGETSENSORVALUE('RP_U3',jSensorLookUpTable,sensorValues)
    ur(1)=VGETSENSORVALUE('RP_UR1',jSensorLookUpTable,sensorValues)
    ur(2)=VGETSENSORVALUE('RP_UR2',jSensorLookUpTable,sensorValues)
    ur(3)=VGETSENSORVALUE('RP_UR3',jSensorLookUpTable,sensorValues)
    dotu=dot_product(u,c)
    s_eff=31.699921242759284d0+dotu-6.0d0*t
    phase_deg=modulo(43200.0d0*t,360.0d0)
    call interpolate_table(s_eff,phase_deg,b,grad,outside)
    if (outside) then
      write(6,*) 'PRECOMPUTED_TABLE_RANGE_ERROR s_eff_mm=',s_eff
      lFlagsDefine(iStopAnalysis)=1
    endif
    ang=sqrt(dot_product(ur,ur))
    if (ang<1.0d-14) then
      moment=m0; a=a0
    else
      kx=ur(1)/ang; ky=ur(2)/ang; kz=ur(3)/ang; cc=cos(ang); ss=sin(ang)
      moment(1)=(cc+kx*kx*(1-cc))*m0(1)+(kx*ky*(1-cc)-kz*ss)*m0(2)+(kx*kz*(1-cc)+ky*ss)*m0(3)
      moment(2)=(ky*kx*(1-cc)+kz*ss)*m0(1)+(cc+ky*ky*(1-cc))*m0(2)+(ky*kz*(1-cc)-kx*ss)*m0(3)
      moment(3)=(kz*kx*(1-cc)-ky*ss)*m0(1)+(kz*ky*(1-cc)+kx*ss)*m0(2)+(cc+kz*kz*(1-cc))*m0(3)
      a(1)=(cc+kx*kx*(1-cc))*a0(1)+(kx*ky*(1-cc)-kz*ss)*a0(2)+(kx*kz*(1-cc)+ky*ss)*a0(3)
      a(2)=(ky*kx*(1-cc)+kz*ss)*a0(1)+(cc+ky*ky*(1-cc))*a0(2)+(ky*kz*(1-cc)-kx*ss)*a0(3)
      a(3)=(kz*kx*(1-cc)-ky*ss)*a0(1)+(kz*ky*(1-cc)+kx*ss)*a0(2)+(cc+kz*kz*(1-cc))*a0(3)
    endif
    b=0.012d0*b; grad=0.002000d0*grad
    if (t<0.0010000d0) then
      q=max(0.0d0,min(1.0d0,t/0.0010000d0)); scale=q*q*(3.0d0-2.0d0*q)
    else
      scale=1.0d0
    endif
    do i=1,3
      loads(i)=scale*sum(moment(:)*grad(:,i))
    enddo
    loads(4)=scale*1000.0d0*(moment(2)*b(3)-moment(3)*b(2))
    loads(5)=scale*1000.0d0*(moment(3)*b(1)-moment(1)*b(3))
    loads(6)=scale*1000.0d0*(moment(1)*b(2)-moment(2)*b(1))

    call open_increment_log()
    write(77,'(9(ES18.10,:,","))') t,s_eff,phase_deg,loads
    last_time=t
  endif
  AmpValueNew=0.0d0
  if (ampName(1:9)=='SOCKET_FX') AmpValueNew=loads(1)
  if (ampName(1:9)=='SOCKET_FY') AmpValueNew=loads(2)
  if (ampName(1:9)=='SOCKET_FZ') AmpValueNew=loads(3)
  if (ampName(1:9)=='SOCKET_MX') AmpValueNew=loads(4)
  if (ampName(1:9)=='SOCKET_MY') AmpValueNew=loads(5)
  if (ampName(1:9)=='SOCKET_MZ') AmpValueNew=loads(6)
  return
end subroutine vuamp

subroutine vexternaldb(lOp,i_Array,niArray,r_Array,nrArray)
  use truecel_precomputed_data
  include 'VABA_PARAM.INC'
  dimension i_Array(niArray),r_Array(nrArray)
  if (lOp==0 .and. .not.table_loaded) call load_magnetic_table()
  if (lOp==6) then
    write(6,*) 'PRECOMPUTED_TABLE_LOOKUPS=',lookup_count,' SOCKET_CALLS=0'
    if (log_open) then
      flush(77); close(77); log_open=.false.
    endif
  endif
  return
end subroutine vexternaldb

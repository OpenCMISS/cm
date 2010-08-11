!> \file
!> $Id:  $
!> \author Sander Land
!> \brief This module contains some hardcoded cell models and integration routines for cardiac electrophysiology.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

module electrophysiology_cell_routines
  USE BASE_ROUTINES
  use field_routines
  use kinds
  use strings
  use types

  implicit none
  private

  !interfaces

  public  bueno_orovio_initialize, bueno_orovio_integrate
  public  tentusscher06_initialize, tentusscher06_integrate

  interface pow
    module procedure powr
    module procedure powi
  end  interface pow

contains
  ! for auto generated code

  real(dp) function powr(a,b)
    real(dp), intent(in) :: a,b
    powr = a**b
  end function powr
  real(dp) function powi(a,b)
    real(dp), intent(in) :: a
    integer(intg), intent(in) :: b
    powi = a**b
  end function powi

  ! Bueno-Orovio (2008) minimal cell model. Membrane potential rescaled, partial evaluation applied.

  ! initialize field
  subroutine bueno_orovio_initialize(field,err,error,*)
    type(field_type), intent(inout), pointer :: field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    integer(intg) :: i
    REAL(DP), PARAMETER, DIMENSION(1:4) :: y0 = (/ -8.09242e+01,  9.99993e-01, 2.16366e-02, 9.84320e-01 /) ! paced initial condition

    call enters('bueno_orovio_init',err,error,*999)
    DO I=1,4
      CALL FIELD_COMPONENT_VALUES_INITIALISE(field,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,y0(I),Err,ERROR,*999)
    END DO

    call exits('bueno_orovio_init')
    return
999 call errors('bueno_orovio_init',err,error)
    call exits('bueno_orovio_init')
    return 1

  end subroutine bueno_orovio_initialize

  subroutine bueno_orovio_dydt(t,y,dydt,activ)
    real(dp), dimension(:), intent(out) :: dydt(:) !< derivatives
    real(dp), dimension(:), intent(in) :: y(:) !< current state
    real(dp), intent(in) :: t  !< current time 
    real(dp), intent(in) :: activ  !< activation factor. 1 for default stimulus current, 0 for none

    real(dp),parameter :: u_m=3.0e-01, u_p=1.3e-01, u_q=6.0e-03, u_r=6.0e-03, vscale=81.0, period = 1000.0
    real(dp) :: J_fi, J_si, J_so, Vm, bv, m, p, q, r, tau_o, tau_s, tau_so, tau_v_minus, tau_w_minus, w_inf, J_stim
  
    if (mod(t,period)>=0 .and. mod(t,period)<=0.999999) then ! apply stimulus in first ms
      J_stim = -0.65 * activ ! stimulus: ~52 mV for 1 ms
    else
     J_stim  = 0
    end if

    bv = (1+(Y(1)/vscale)); !< rescale, original Y(1) in model
    if (bv < u_m) then
      m = 0
    else
      m = 1
    end if
    if (bv < u_q) then
      q = 0
    else
      q = 1
    end if
    if (bv < u_p) then
      p = 0
    else
      p = 1
    end if

    if (bv < u_r) then
      r = 0
    else
      r = 1
    end if

    J_fi = (9.09090909090909e+00*(-m)*Y(2)*(-3.0e-01+bv)*(1.55e+00-bv))
    tau_v_minus = ((1150*q)+(60*(1-q)))
    dydt(2) = (((1-m)*(0-Y(2))/tau_v_minus)-(6.89655172413793e-01*m*Y(2)))
    tau_o = ((400*(1-r))+(6*r))
    tau_so = (3.002e+01+(-1.4512e+01*(1+tanh((2.046e+00*(-6.5e-01+bv))))))
    J_so = ((bv*(1-p)/tau_o)+(p/tau_so))
    J_si = (5.29801324503311e-01*(-p)*Y(4)*Y(3))
    dydt(1) = ((-J_fi-J_so-J_si-J_stim)*vscale)
    Vm = (-83+(8.57e+01*bv))
    tau_s = ((2.7342e+00*(1-p))+(16*p))
    dydt(3) = (((5.0e-01*(1+tanh((2.0994e+00*(-9.087e-01+bv)))))-Y(3))/tau_s)
    w_inf = (((1-r)*(1-(1.42857142857143e+01*bv)))+(9.4e-01*r))
    tau_w_minus = (60+(-2.25e+01*(1+tanh((65*(-3.0e-02+bv))))))
    dydt(4) = (((1-r)*(w_inf-Y(4))/tau_w_minus)-(5.0e-03*r*Y(4)))
  end subroutine bueno_orovio_dydt   

  ! integrates all cell models
  subroutine bueno_orovio_integrate(cells,materials,t0,t1,err,error,*)
    type(field_type), intent(inout), pointer :: cells, materials !<Independent field storing the cell data, and material field component 1 the activation flags
    real(dp), intent(in)    :: t0, t1 !<Integrate from time t0 to t1 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    integer, parameter :: celldim = 4
    real(dp), dimension(1:celldim) :: y, dydt
    real(dp) :: t, dt, activ
    real(dp), dimension(:), pointer :: celldata, activdata

    integer(intg) :: ncells, i, d, nodeno
    type(domain_ptr_type), pointer :: domain
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: CELLS_VARIABLE, ACTIV_VARIABLE

    call enters('bueno_orovio_integrate',err,error,*999)

    domain=>cells%decomposition%domain(1)
    ncells = domain%ptr%topology%nodes%number_of_nodes ! local

    CELLS_VARIABLE=>cells%VARIABLE_TYPE_MAP(field_u_variable_type)%PTR
    ACTIV_VARIABLE=>materials%VARIABLE_TYPE_MAP(field_u_variable_type)%PTR

    CALL FIELD_PARAMETER_SET_DATA_GET(cells,field_u_variable_type,field_values_set_type,celldata,ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SET_DATA_GET(materials,field_u_variable_type,field_values_set_type,activdata,ERR,ERROR,*999)

    do i=1,ncells
      !   field ->   y
      do d=1,celldim
        y(d) = celldata(CELLS_VARIABLE%COMPONENTS(d)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,i))
      end do
      activ = activdata(ACTIV_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,i))
      
      t = t0
      do while (t < t1 - 1e-6)
        ! integrate one cell, one time step. 
        ! just use adaptive forward euler, tested in matlab
        call bueno_orovio_dydt(t,y,dydt,activ)
        dt = min(t1-t,0.5)  ! default max
        dt = min(dt,1 / abs(dydt(1))) ! maximum increase in Vm
!        IF(DT < 1e-6)  WRITE(*,*) 'CELL ',i, 'DVDT = ',dydt(1),'dt = ',dt,' activ=',activ,'y=',y
        y = y + dt * dydt ! FWE
        t = t + dt
      end do
      !   y -> field  
      do d=1,celldim
 !       call field_parameter_set_update_local_node(cells,field_u_variable_type,field_values_set_type,1,i,d,y(d), err,error,*999)
        celldata(CELLS_VARIABLE%COMPONENTS(d)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,i)) = y(d)
      end do
    end do
    CALL FIELD_PARAMETER_SET_DATA_RESTORE(cells,field_u_variable_type,field_values_set_type,celldata,ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SET_DATA_RESTORE(materials,field_u_variable_type,field_values_set_type,activdata,ERR,ERROR,*999)


    call exits('bueno_orovio_integrate')
    return
999 call errors('bueno_orovio_integrate',err,error)
    call exits('bueno_orovio_integrate')
    return 1
  end subroutine bueno_orovio_integrate
  


  ! Ten Tusscher and Panfilov (2006) human ventricular cell model

  ! initialize field
  subroutine tentusscher06_initialize(field,err,error,*)
    type(field_type), intent(inout), pointer :: field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    integer(intg) :: i
    REAL(DP), PARAMETER, DIMENSION(1:19) :: y0 = (/ -85.23, 0.9755, 0.9953, 0.7888, 3.64, 0.000126, 0.00036, 0.9073, 0.7444,&
        & 0.7045, 0.00172,  3.373e-5, 136.89, 0.00621, 0.4712, 0.0095, 8.604, 2.42e-8, 0.999998 /) ! paced initial condition

    call enters('tentusscher06_init',err,error,*999)
    DO I=1,19
      CALL FIELD_COMPONENT_VALUES_INITIALISE(field,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,y0(I),Err,ERROR,*999)
    END DO

    call exits('tentusscher06_init')
    return
999 call errors('tentusscher06_init',err,error)
    call exits('tentusscher06_init')
    return 1

  end subroutine tentusscher06_initialize

  subroutine tentusscher06_dydt(t,y,dydt,activ)
    real(dp), dimension(:), intent(out) :: dydt(:) !< derivatives
    real(dp), dimension(:), intent(in) :: y(:) !< current state
    real(dp), intent(in) :: t  !< current time 
    real(dp), intent(in) :: activ  !< activation factor. 1 for default stimulus current, 0 for none

    real(dp),parameter :: period=1000, activate=1
    real(dp) :: Ca_i_bufc, Ca_sr_bufsr, Ca_ss_bufss, E_Ca, E_K, E_Ks, E_Na, O, alpha_K1, alpha_d, alpha_m, alpha_xr1, alpha_xr2,&
    & alpha_xs, beta_K1, beta_d, beta_m, beta_xr1, beta_xr2, beta_xs, d_inf, f2_inf, fCass_inf, f_inf, gamma_d, h_inf, i_CaL,&
    & i_K1, i_Kr, i_Ks, i_Na, i_NaCa, i_NaK, i_Stim, i_b_Ca, i_b_Na, i_leak, i_p_Ca, i_p_K, i_rel, i_to, i_up, i_xfer, j_inf, k1,&
    & k2, kcasr, m_inf, r_inf, s_inf, tau_d, tau_f, tau_f2, tau_fCass, tau_h, tau_j, tau_m, tau_r, tau_s, tau_xr1, tau_xr2,tau_xs,&
    & xK1_inf, xr1_inf, xr2_inf, xs_inf

    if (mod(t,period)>=0 .and. mod(t,period)<=1.999999) then ! apply stimulus in 2 ms
      i_Stim = -3.57142857142857e+01 * activ
    else
      i_Stim = 0
    end if

    i_CaL = (5.75002020961245e-01*Y(12)*Y(4)*Y(2)*Y(3)*(-15+Y(1))*(-2+(2.5e-01*Y(7)*exp((7.48677816454909e-02*(-15+Y(1))))))/&
          & (-1+exp((7.48677816454909e-02*(-15+Y(1))))))
    
    d_inf = (1/(1+exp((1.33333333333333e-01*(-8-Y(1))))))
    alpha_d = (2.5e-01+(1.4e+00/(1+exp((7.69230769230769e-02*(-35-Y(1)))))))
    beta_d = (1.4e+00/(1+exp((2.0e-01*(5+Y(1))))))
    gamma_d = (1/(1+exp((5.0e-02*(50-Y(1))))))
    tau_d = ((1*alpha_d*beta_d)+gamma_d)
    f2_inf = (3.3e-01+(6.7e-01/(1+exp((1.42857142857143e-01*(35+Y(1)))))))
    tau_f2 = ((562*exp((4.16666666666667e-03*(-pow((27+Y(1)),2)))))+(31/(1+exp((1.0e-01*(25-Y(1))))))+&
           & (80/(1+exp((1.0e-01*(30+Y(1)))))))
    dydt(2) = ((f2_inf-Y(2))/tau_f2)
    fCass_inf = (4.0e-01+(6.0e-01/(1+pow((20*Y(7)),2))))
    tau_fCass = (2+(80/(1+pow((20*Y(7)),2))))
    dydt(3) = ((fCass_inf-Y(3))/tau_fCass)
    f_inf = (1/(1+exp((1.42857142857143e-01*(20+Y(1))))))
    tau_f = (20+(1.1025e+03*exp((4.44444444444444e-03*(-pow((27+Y(1)),2)))))+(200/(1+exp((1.0e-01*(13-Y(1))))))+&
          &(180/(1+exp((1.0e-01*(30+Y(1)))))))
    dydt(4) = ((f_inf-Y(4))/tau_f)
    E_Ca = (1.33568803298478e+01*log((2/Y(6))))
    i_b_Ca = (5.92e-04*(Y(1)-E_Ca))
    kcasr = (2.5e+00-(1.5e+00/(1+pow((1.5e+00/Y(5)),2))))
    k1 = (1.5e-01/kcasr)
    O = (k1*pow(Y(7),2)*Y(8)/(6.0e-02+(k1*pow(Y(7),2))))
    i_rel = (1.02e-01*O*(Y(5)-Y(7)))
    i_up = (6.375e-03/(1+(6.25e-08/pow(Y(6),2))))
    i_leak = (3.6e-04*(Y(5)-Y(6)))
    i_xfer = (3.8e-03*(Y(7)-Y(6)))
    k2 = (4.5e-02*kcasr)
    dydt(8) = (((-k2)*Y(7)*Y(8))+(5.0e-03*(1-Y(8))))
    Ca_i_bufc = (1/(1+(2.0e-04/pow((1.0e-03+Y(6)),2))))
    Ca_sr_bufsr = (1/(1+(3/pow((3.0e-01+Y(5)),2))))
    Ca_ss_bufss = (1/(1+(1.0e-04/pow((2.5e-04+Y(7)),2))))
    i_p_Ca = (1.238e-01*Y(6)/(5.0e-04+Y(6)))
    i_NaCa = (8.66622022994244e-05*((2*exp((1.31018617879609e-02*Y(1)))*pow(Y(17),3))-(6860000*&
           &exp((-2.43320290347846e-02*Y(1)))*Y(6)))/(1+(1.0e-01*exp((-2.43320290347846e-02*Y(1))))))
    dydt(6) = (Ca_i_bufc*((6.66910509631797e-02*(i_leak-i_up))+i_xfer-(5.84427487220097e-05*(i_b_Ca+i_p_Ca-(2*i_NaCa)))))
    dydt(5) = (Ca_sr_bufsr*(i_up-i_rel-i_leak))
    dydt(7) = (Ca_ss_bufss*((-1.75328246166029e-02*i_CaL)+(2.00073152889539e+01*i_rel)-(300*i_xfer)))
    E_Na = (2.67137606596956e+01*log((140/Y(17))))
    i_Na = (1.4838e+01*pow(Y(11),3)*Y(9)*Y(10)*(Y(1)-E_Na))
    h_inf = (1/pow((1+exp((1.34589502018843e-01*(7.155e+01+Y(1))))),2))
    if (Y(1) < -40.0) then
      tau_h = (1/((5.7e-02*exp((1.47058823529412e-01*(-80-Y(1)))))+(2.7e+00*exp((7.9e-02*Y(1))))+(310000*exp((3.485e-01*Y(1))))))
    else
      tau_h = (1.68831168831169e-01*(1+exp((-9.00900900900901e-02*(1.066e+01+Y(1))))))
    end if
    dydt(9) = ((h_inf-Y(9))/tau_h)
    j_inf = (1/pow((1+exp((1.34589502018843e-01*(7.155e+01+Y(1))))),2))
    if (Y(1) < -40.0) then
     tau_j = (1/((1*((-25428*exp((2.444e-01*Y(1))))-(6.948e-06*exp((-4.391e-02*Y(1)))))*(3.778e+01+Y(1))/&
          &(1+exp((3.11e-01*(7.923e+01+Y(1))))))+(2.424e-02*exp((-1.052e-02*Y(1)))/(1+exp((-1.378e-01*(4.014e+01+Y(1))))))))
    else
     tau_j = (1.66666666666667e+00/exp((5.7e-02*Y(1)))*(1+exp((-1.0e-01*(32+Y(1))))))
    end if
    dydt(10) = ((j_inf-Y(10))/tau_j)
    m_inf = (1/pow((1+exp((1.10741971207087e-01*(-5.686e+01-Y(1))))),2))
    alpha_m = (1/(1+exp((2.0e-01*(-60-Y(1))))))
    beta_m = ((1.0e-01/(1+exp((2.0e-01*(35+Y(1))))))+(1.0e-01/(1+exp((5.0e-03*(-50+Y(1)))))))
    tau_m = (1*alpha_m*beta_m)
    dydt(11) = ((m_inf-Y(11))/tau_m)
    E_K = (2.67137606596956e+01*log((5.4e+00/Y(13))))
    alpha_K1 = (1.0e-01/(1+exp((6.0e-02*(-200+Y(1)-E_K)))))
    beta_K1 = (((3*exp((2.0e-04*(100+Y(1)-E_K))))+exp((1.0e-01*(-10+Y(1)-E_K))))/(1+exp((-5.0e-01*(Y(1)-E_K)))))
    xK1_inf = (alpha_K1/(alpha_K1+beta_K1))
    i_K1 = (5.405e+00*xK1_inf*(Y(1)-E_K))
    i_to = (2.94e-01*Y(18)*Y(19)*(Y(1)-E_K))
    i_Kr = (1.53e-01*Y(14)*Y(15)*(Y(1)-E_K))
    E_Ks = (2.67137606596956e+01*log((9.6e+00/(Y(13)+(3.0e-02*Y(17))))))
    i_Ks = (3.92e-01*pow(Y(16),2)*(Y(1)-E_Ks))
    i_NaK = (2.298375e+00*Y(17)/(40+Y(17))/(1+(1.245e-01*exp((-3.74338908227455e-03*Y(1))))+&
          &(3.53e-02*exp((3.74338908227455e-02*(-Y(1)))))))
    i_b_Na = (2.9e-04*(Y(1)-E_Na))
    i_p_K = (1.46e-02*(Y(1)-E_K)/(1+exp((1.67224080267559e-01*(25-Y(1))))))
    dydt(12) = ((d_inf-Y(12))/tau_d)
    dydt(1) = (-1*(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_b_Na+i_NaCa+i_b_Ca+i_p_K+i_p_Ca+i_Stim))
    dydt(13) = (-1.16885497444019e-04*(i_K1+i_to+i_Kr+i_Ks+i_p_K+i_Stim-(2*i_NaK)))
    xr1_inf = (1/(1+exp((1.42857142857143e-01*(-26-Y(1))))))
    alpha_xr1 = (450/(1+exp((1.0e-01*(-45-Y(1))))))
    beta_xr1 = (6/(1+exp((8.69565217391304e-02*(30+Y(1))))))
    tau_xr1 = (1*alpha_xr1*beta_xr1)
    dydt(14) = ((xr1_inf-Y(14))/tau_xr1)
    xr2_inf = (1/(1+exp((4.16666666666667e-02*(88+Y(1))))))
    alpha_xr2 = (3/(1+exp((5.0e-02*(-60-Y(1))))))
    beta_xr2 = (1.12e+00/(1+exp((5.0e-02*(-60+Y(1))))))
    tau_xr2 = (1*alpha_xr2*beta_xr2)
    dydt(15) = ((xr2_inf-Y(15))/tau_xr2)
    xs_inf = (1/(1+exp((7.14285714285714e-02*(-5-Y(1))))))
    alpha_xs = (1400/sqrt((1+exp((1.66666666666667e-01*(5-Y(1)))))))
    beta_xs = (1/(1+exp((6.66666666666667e-02*(-35+Y(1))))))
    tau_xs = (80+(1*alpha_xs*beta_xs))
    dydt(16) = ((xs_inf-Y(16))/tau_xs)
    dydt(17) = (-1.16885497444019e-04*(i_Na+i_b_Na+(3*i_NaK)+(3*i_NaCa)))
    r_inf = (1/(1+exp((1.66666666666667e-01*(20-Y(1))))))
    tau_r = (8.0e-01+(9.5e+00*exp((5.55555555555556e-04*(-pow((40+Y(1)),2))))))
    dydt(18) = ((r_inf-Y(18))/tau_r)
    s_inf = (1/(1+exp((2.0e-01*(20+Y(1))))))
    tau_s = (3+(85*exp((3.125e-03*(-pow((45+Y(1)),2)))))+(5/(1+exp((2.0e-01*(-20+Y(1)))))))
    dydt(19) = ((s_inf-Y(19))/tau_s)
  end subroutine tentusscher06_dydt   

  ! integrates all cell models
  subroutine tentusscher06_integrate(cells,materials,t0,t1,err,error,*)
    type(field_type), intent(inout), pointer :: cells, materials !<Independent field storing the cell data, and material field component 1 the activation flags
    real(dp), intent(in)    :: t0, t1 !<Integrate from time t0 to t1 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    integer, parameter :: celldim = 19
    real(dp), dimension(1:celldim) :: dydt
    real(dp), dimension(:), pointer :: y
    real(dp) :: t, dt, activ, m_inf, d_inf, m_inf0, d_inf0
    real(dp), dimension(:), pointer :: celldata, activdata

    integer(intg) :: ncells, i, d, nodeno
    type(domain_ptr_type), pointer :: domain
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: CELLS_VARIABLE, ACTIV_VARIABLE

    call enters('tentusscher06_integrate',err,error,*999)

    domain=>cells%decomposition%domain(1)
    ncells = domain%ptr%topology%nodes%number_of_nodes ! local

    CELLS_VARIABLE=>cells%VARIABLE_TYPE_MAP(field_u_variable_type)%PTR
    ACTIV_VARIABLE=>materials%VARIABLE_TYPE_MAP(field_u_variable_type)%PTR

    CALL FIELD_PARAMETER_SET_DATA_GET(cells,field_u_variable_type,field_values_set_type,celldata,ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SET_DATA_GET(materials,field_u_variable_type,field_values_set_type,activdata,ERR,ERROR,*999)

    do i=1,ncells
      d = CELLS_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,i)
      y => celldata(d:d+celldim-1)
      !   field ->   y
  !    do d=1,celldim
!        nodeno = domain%ptr%topology%nodes%nodes(i)%global_number
!        call field_parameter_set_get_node(cells,field_u_variable_type,field_values_set_type,1,nodeno,d,y(d),err,error,*999)
 !       y(d) = celldata(CELLS_VARIABLE%COMPONENTS(d)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,i))
 !     end do
      ! field_parameter_set_get_node(materials,field_u_variable_type,field_values_set_type,1,nodeno,1,activ,err,error,*999)
      activ = activdata(ACTIV_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,i))
      
  
      t = t0
      do while (t < t1 - 1e-6)
        ! integrate one cell, one time step. 
        ! just use adaptive forward euler, tested in matlab
        call tentusscher06_dydt(t,y,dydt,activ)
        dt = min(t1-t,1.0/3)  ! default max
        dt = min(dt,1 / abs(dydt(1))) ! maximum increase in Vm
        
        m_inf0 = (1/pow((1+exp((1.10741971207087e-01*(-5.686e+01-Y(1))))),2))
        d_inf0 = (1/(1+exp((1.33333333333333e-01*(-8-Y(1))))))
        
        y = y + dt * dydt ! FWE
        
        m_inf = (1/pow((1+exp((1.10741971207087e-01*(-5.686e+01-Y(1))))),2))
        d_inf = (1/(1+exp((1.33333333333333e-01*(-8-Y(1))))))
        
        ! FWEoo on m,d gates : prevent gates from crossing midpoint steady state
        
        m_inf = (m_inf0 + m_inf)/2
        d_inf = (d_inf0 + d_inf)/2
        if(dydt(11) * (y(11) - m_inf) > 0) then
          y(11) = m_inf
        end if
        if(dydt(12) * (y(12) - d_inf) > 0) then
          y(12) = d_inf
        end if        
        
        t = t + dt
      end do
      !   y -> field  
      !do d=1,celldim
 !       call field_parameter_set_update_local_node(cells,field_u_variable_type,field_values_set_type,1,i,d,y(d), err,error,*999)
     !   celldata(CELLS_VARIABLE%COMPONENTS(d)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,i)) = y(d)
    !  end do
    end do
    CALL FIELD_PARAMETER_SET_DATA_RESTORE(cells,field_u_variable_type,field_values_set_type,celldata,ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SET_DATA_RESTORE(materials,field_u_variable_type,field_values_set_type,activdata,ERR,ERROR,*999)

    call exits('tentusscher06_integrate')
    return
999 call errors('tentusscher06_integrate',err,error)
    call exits('tentusscher06_integrate')
    return 1
  end subroutine tentusscher06_integrate
end module electrophysiology_cell_routines



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

contains
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

  ! Bueno-Orovio (2008) minimal cell model. Membrane potential rescaled, partial evaluation applied.
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
    type(field_type), intent(inout), pointer :: cells, materials
    real(dp), intent(in)    :: t0, t1
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    integer, parameter :: celldim = 4
    real(dp), dimension(1:celldim) :: y, dydt
    real(dp) :: t, dt, activ

    integer(intg) :: ncells, i, d, nodeno
    type(domain_ptr_type), pointer :: domain

    call enters('bueno_orovio_integrate',err,error,*999)

    domain=>cells%decomposition%domain(1)  ! ?
    ncells = domain%ptr%topology%nodes%number_of_nodes ! local

    do i=1,ncells
      !   field ->   y
      do d=1,celldim
        nodeno = domain%ptr%topology%nodes%nodes(i)%global_number
        call field_parameter_set_get_node(cells,field_u_variable_type,field_values_set_type,1,nodeno,d,y(d),err,error,*999)
      end do
      call field_parameter_set_get_node(materials,field_u_variable_type,field_values_set_type,1,nodeno,1,activ,err,error,*999)

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
        call field_parameter_set_update_local_node(cells,field_u_variable_type,field_values_set_type,1,i,d,y(d), err,error,*999)
      end do
    end do

    call exits('bueno_orovio_integrate')
    return
999 call errors('bueno_orovio_integrate',err,error)
    call exits('bueno_orovio_integrate')
    return 1
  end subroutine bueno_orovio_integrate
  
end module electrophysiology_cell_routines



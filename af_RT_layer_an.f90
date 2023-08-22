!=======================================================================================================
! The subroutine calculates probability of reflection, penetration through and absorption in
! a layer of given optical thickness tau_0.
! Result is approximate and corresponds to tau>>1!
!   lambda - probability for photon to survive after a scattering: lambda=1. for conservative scattering.
!=======================================================================================================
subroutine RT_prob_layer_asym(tau_0,lambda,P_ref,P_abs,P_pen)
implicit none
real*8,intent(in)::tau_0,lambda
real*8::k,tau_e,P_ref,P_pen,P_abs
  k=get_k(lambda,1.d-3)
  tau_e=get_tau_e(lambda)
  !==the basic propabilities=============================!
  P_ref=1.d0 - sqrt(1.d0-lambda)/tanh(k*(tau_0+2*tau_e))
  P_pen=sqrt(1.d0-lambda)/sinh(k*(tau_0+2*tau_e))
  P_abs=1.d0-P_pen-P_ref
  !======================================================!
return
contains
  !===================================================================
  ! This function is used to get the propabilities in RT of a layer.
  ! DIN (59)
  ! works badly for small lambda (lambda<0.5)
  !===================================================================
  real*8 function get_tau_e(lambda)
  implicit none
  real*8,intent(in)::lambda
  real*8::k,a,b
  real*8::pi=3.141592653589793d0
  real*8::int_simpson_1of3   !==functions==!
    k=get_k(lambda,1.d-3)
    a=0.d0 !1.d-9
    b=1.d0-a
    get_tau_e=int_simpson_1of3(fun,a,1.d-2,1.d-2,lambda,k)&
             +int_simpson_1of3(fun,1.d-2,0.99d0,1.d-2,lambda,k)&
             +int_simpson_1of3(fun,0.99d0,b,1.d-2,lambda,k)
    get_tau_e=1.d0/lambda - get_tau_e/pi
  return
  end function get_tau_e
  !===================================================================

  real*8 function fun(nu,lambda,k)
  implicit none
  real*8,intent(in)::nu,lambda,k
  real*8::pi=3.141592653589793d0
  real*8::Xi
    !Xi=nu/2*log( (1.d0+nu)/(1.d0-nu) )
    Xi=( 1.d0/lambda-nu/2*log( (1.d0+nu)/(1.d0-nu) ) )/pi/nu*2
    fun=atan(Xi)/(1.d0-k**2*nu**2)
  return
  end function fun

  !=================================================================
  ! This function is used to get the propabilities in RT of a layer.
  ! DIN (56).
  ! function gets characteristic number .
  ! lambda\in (0.,1.)
  !=================================================================
  real*8 function get_k(lambda,eps)
  implicit none
  real*8,intent(in)::lambda,eps
  real*8::lambda_,res,k1,k2,k
  integer::det
    lambda_=1.d0/lambda
    det=1
    k1=0.d0; k2=0.d0
    do while(det.eq.1)
      k2=k1+(1.d0-k1)/2
      if(rhs(k2).ge.lambda_)then
        det=2
      else
        k1=k2
      end if
    end do
    k=k1
    do while((k2-k1).ge.eps)
      k=(k1+k2)/2
      if(rhs(k).le.lambda_)then
        k1=k
      else
        k2=k
      end if
    end do
    get_k=k
  return
  end function get_k

  real*8 function rhs(k)
  implicit none
  real*8::k
    rhs=0.5d0/k*log((1.d0+k)/(1.d0-k))
  return
  end function rhs
end subroutine RT_prob_layer_asym


subroutine test_RT_prob_layer_asym()
implicit none
real*8::lambda,tau_0,P_pen,P_abs,P_ref
  tau_0=0.001d0
  lambda=0.9d0       !== probability of photon to survive in a single scattering event ==!
  call RT_prob_layer_asym(tau_0,lambda,P_ref,P_abs,P_pen)
  write(*,*)tau_0,lambda,P_ref,P_abs,P_pen
  lambda=0.9999d0
  call RT_prob_layer_asym(tau_0,lambda,P_ref,P_abs,P_pen)
  write(*,*)tau_0,lambda,P_ref,P_abs,P_pen
return
end subroutine test_RT_prob_layer_asym
!=================================================================================================



!=================================================================================================
! The subroutine calculates the array containing data on probability of reflection from the layer.
!   lambda - probability of photon to survive in a single scattering event.
!   [tau_min,tau_max] - the considered interval of optical thickness.
!=================================================================================================
subroutine get_array_P_ref(mas_P_ref,n,tau_min,tau_max,lambda)
implicit none
real*8::mas_P_ref
dimension mas_P_ref(n,2)
integer,intent(in)::n
real*8,intent(in)::tau_min,tau_max,lambda
integer::i
real*8::tau,dtau,P_ref,P_abs,P_pen
  dtau=(tau_max-tau_min)/(n-1)
  i=1
  do while(i.le.n)
    tau=tau_min+(i-1)*dtau
    call RT_prob_layer_asym(tau,lambda,P_ref,P_abs,P_pen)
    mas_P_ref(i,1)=tau
    mas_P_ref(i,2)=P_ref
    !write(*,*)tau,P_ref,P_abs,P_pen; read(*,*)
    i=i+1
  end do
return
end subroutine get_array_P_ref
!=================================================================================================

!=======================================================================
!=======================================================================
module mod_sig_ave
  integer::l_pol
  real*8::sig1,sig2,bb,TT,mas_,knt
  dimension sig1(10,300,2), sig2(10,300,2), mas_(300,2),knt(100,2)
end module mod_sig_ave


!======================================================================
!  Test routine for maxL: calculates luminosities for C, B, and A zones
!  as a function of magnetic field parameter b and column height H.
!======================================================================
subroutine test_maxL_CBA()
implicit none
real*8 :: b, H
real*8 :: resC, resB, resA
real*8 :: z2RmC, z2RmB, z2RmA
integer :: iC, iB, iA
integer :: jC, jB, jA

  b = 2.d0
  write(*,*) "# CBA-zones"
  write(*,*) "# length = 1"

  do while (b .le. 22.2d0)
    write(*,*) "# b=", b
    H = 0.1d0
    do while (H .le. 1.5d0)
      call maxL(b, H, 1.4d0, 10.d0, resC, iC, jC, z2RmC, 1)
      !call maxL(b, H, 1.4d0, 10.d0, resB, iB, jB, z2RmB, 2)
      !call maxL(b, H, 1.4d0, 10.d0, resA, iA, jA, z2RmA, 3)
      write(*,*) b, H, resC, z2RmC!, resB, z2RmB, resA, z2RmA
      H = H * 1.5d0
    end do
    b = b * 10.d0
  end do
return
end subroutine test_maxL_CBA


!=================================================================================================================
!  The subroutine calculates the luminosity which corresponds to the accretion column
!  of height H. [H]=[h/R]
!  This is 2-d code, where the accretion column structure across the field lines is taken into account.
!  b0 - surface magnetic field strength in units of critical B-field strength
!  M - neutron star mass measured in solar masses, R - [NS radius]/1[km]
!  det_S defines the zone, where the accretion disc is interrupted: 1 - C-zone, 2 - B-zone, 3 - A-zone
!  Output: res=L39, det_T and det_nu check if the temperature is high enough to break the BF and create e-e+ pairs
!          z2Rm - z/Rm
!=================================================================================================================
subroutine maxL(b0,H,M,R,res,det_T,det_nu,z2Rm,det_S)
implicit none
real*8,intent(in)::b0,H,M,R
integer,intent(in)::det_S
real*8,intent(out)::res,z2Rm
integer,intent(out)::det_T,det_nu

real*8::res1,rel_err
integer::n,det

  det=13
  n=10
  call maxL_n(b0,H,M,R,n,res,det_T,det_nu,z2Rm,det_S)

  do while((det.eq.13).and.(n.le.32))
    n=n*2
    call maxL_n(b0,H,M,R,n,res1,det_T,det_nu,z2Rm,det_S)
    if(abs(res1).gt.0.d0)then
      rel_err=abs((res-res1)/res1)
    else
      rel_err=abs(res-res1)
    end if
    if(rel_err.le.0.03d0)then
      det=1
    end if
    res=res1
    !write(*,*)n,res,rel_err
  end do
return
end subroutine maxL


!============================================================================================
!  Calculates the luminosity L39 corresponding to a given accretion column height H=h/R
!  for a fixed grid resolution nn and a selected disc-interruption zone det_S.
!  The routine iterates over L39, updates the column geometry, temperature, opacity,
!  and integrates the emergent flux over the column height.
!============================================================================================
subroutine maxL_n(b0,H,M,R,nn,res,det_T,det_nu,z2Rm,det_S)
implicit none
real*8,intent(in)::b0,H,M,R
integer,intent(in)::nn
integer,intent(in)::det_S !==what zone of the accretion disc which we use==!
real*8,intent(out)::res,z2Rm
integer,intent(out)::det_T,det_nu
integer::i,j,itt
real*8::hh,dh,res1,Rm,z,alpha_disc,P_g_rad,Teff
real*8::masRes,T_new,L39,B12,Lambda,T_corr,l2d,length,S_d,d,alpha,T_ort,b_h
real*8::vel_top,I1_h,I2_0
real*8::columnVel,sigma_average_ross_ort,sigma_average_ross_par,I1,I2,find_base   !==functions==!
real*8::pi
dimension masRes(2*nn+1),T_new(2*nn+1),T_ort(10,2)

  pi=3.14159265359d0

  alpha=1.d0    !==velocity coefficient: power which describes the velocity evolution with the height==!
  alpha_disc=0.1d0
  Lambda=0.5d0
  T_corr=1.d0
  length=1.d0 !0.5d0

  B12=b0*44.13d0

  det_nu=0
  det_T=0

  T_new(1:2*nn+1)=30.d0

  res=0.5d0
  L39=0.1d0
  hh=0.001d0
!write(*,*)"#a1"

  itt=1    !==iterative step==!

  do while((abs((res-L39)/L39).gt.0.02d0).and.(itt.le.25))
    !write(*,*)"#a2",abs((res-L39)/L39),L39,H
    L39=(res+L39)/2.d0

111 Rm=7.d7*Lambda*(M**(1.d0/7.d0))*(B12**(4.d0/7.d0))*(L39**(-2.d0/7.d0)) !; write(*,*)"*** ",L39,Rm
    res=0.d0
    res1=0.d0
    dh=(H-hh)/(2.d0*nn)  !dh=(H)/(2.d0*nn)
    masRes(1:2*nn+1)=0.d0

    !==define accretion channel geometry==!
    select case(det_S)
      case(1)
        z=Rm*0.056d0*(L39**(3.d0/20.d0))/(M**(21.d0/40.d0))*((Rm/1.d8)**(1.d0/8.d0))*(alpha_disc**(-0.1d0))     !==C-zone==!
      case(2)
        z=Rm*0.064d0*(L39**(1.d0/5.d0))/(M**(11.d0/20.d0))*((Rm/1.d8)**(1.d0/20.d0))*(alpha_disc**(-0.1d0))     !==B-zone==!
      case(3)
        z=1.d7*L39/M                                                                                            !==A-zone==!
      case default
        z=Rm*0.056d0*(L39**(3.d0/20.d0))/(M**(21.d0/40.d0))*((Rm/1.d8)**(1.d0/8.d0))*(alpha_disc**(-0.1d0))     !==C-zone==!
    end select

    z=z/2.d0

    d=1.d6*abs(asin(sqrt(1.d6/Rm))-asin(sqrt(1.d6/(Rm-2.d0*z))))
    l2d=2.d0*pi*length*1.d6*sqrt(1.d6/(Rm-z))/d
    S_d=2.d0*pi*length*1.d6*sqrt(1.d6/(Rm-z))*d

    !d=1.3d4; l=7.7d5; l2d=l/d; S_d=l*d
    !=====================================!

    if(((Rm-2.d0*z).le.0.d0).or.(abs(1.d6/(Rm-2.d0*z)).ge.1.d0))then
      L39=L39/1.1d0
      goto 111
    end if

    vel_top=columnVel(H*R,H*R,M,R,alpha)

    i=1

    do while(i.le.(2*nn+1))
      b_h=b0*((R/(R+R*hh))**3)
      I1_h=I1(hh,H,alpha)

      !==first estimate of masRes; mostly diagnostic/legacy, overwritten below==!
      masRes(i)=I1_h
      !masRes(i)=masRes(i)*columnVel(hh*R,H*R,M,R,alpha)/vel_top*(((R+hh)/R)**3)/&
      !          sigma_average_ross_ort(b_h,T_new(i)*T_corr)!*1.17d0*1.17d0   !==1.17 because of chemical composition==!
masRes(i)=masRes(i)*columnVel(hh*R,H*R,M,R,alpha)/vel_top*(((R+hh)/R)**3)/&
          1.d0 !*1.17d0*1.17d0   !==1.17 because of chemical composition==!

      !T_new(i)=(6.64d14*L39/S_d/vel_top*abs(I1_h))**0.25d0
      !T_new(i)=(6.64d14*L39/S_d/vel_top*abs(I1_h)+&
      !          9.d0*M/sigma_average_ross_par(b0*((R/(R+R*H))**3),T_new(2*nn+1))/((R+R*H)/10.d0)**2)**0.25d0
      T_new(i)=(6.64d14*L39/S_d/vel_top*abs(I1_h)+&
                9.d0*M/1.d0/((R+R*H)/10.d0)**2)**0.25d0

      if(T_new(i).le.0.d0)then
        T_new(i)=1.d-4
      end if

      j=1
      do while(j.le.10)
        T_ort(j,1)=(j-1)*1.d0/9.d0
        j=j+1
      end do

      T_ort(1:10,2)=T_new(i)

      !==calculations of the accretion column structure across the magnetic field are switched off==!
      j=2
      goto 333

      !==calculations of the accretion column structure across the magnetic field==!
      do while(j.le.10)
        T_ort(j,2)=(6.64d14*L39/S_d/vel_top*abs(I1_h)*&
                    I2(b_h,(j-1)*1.d0/9.d0,T_ort,10)/I2(b_h,0.d0,T_ort,10))**0.25d0
        !write(*,*)j,T_ort(1,2),T_ort(j,2)
        j=j+1
      end do

333   j=j

      I2_0=I2(b_h,0.d0,T_ort,10)

      masRes(i)=I1_h/I2_0
      masRes(i)=masRes(i)*columnVel(hh*R,H*R,M,R,alpha)/vel_top*(((R+hh)/R)**3)

      !===temperature checking=====================
      if(T_new(i).gt.(1300.d0*sqrt(b_h)))then
        !masRes(i)=0.d0
        if(T_new(i).ge.1.d3)then
          det_nu=2
        end if
      else
        if(T_new(i).ge.1.d3)then
          det_nu=1
        end if
      end if

      P_g_rad=3.d20*(L39**0.6d0)*(B12**(-0.5d0))/columnVel(hh*R,H*R,M,R,alpha)*1.5d-9*T_new(i)&
              +5.2d13/3.d0*T_new(i)**4

      if(P_g_rad/b_h/b_h/1.d24*4.d0*pi.ge.1.d0)then
        det_T=1
      end if
      !============================================

      res1=res1+dh*4.d0*M*(l2d/50.d0)*masRes(i)

      Teff=((1.d39/25.d0*M/1.d6/d*masRes(i)/5.67d-5)**(0.25d0))/11604.d0/1.d3/1.22d0

      !!write(*,*)masRes(i),columnVel(hh*R,H*R,M,R,alpha)
      !write(*,*)i,b0,H,hh,T_new(i)!,sigma_average_ross_ort(b_h,T_new(i)*T_corr),Teff,columnVel(hh*R,H*R,M,R,alpha)!, det_T !P_g_rad/b_h/1.d24*4.d0*pi !det_T!  ,l2d        !,det_S
      i=i+1
      hh=hh+dh
    end do

    i=1
    do while(i.le.(2*nn-1))
      res=res+2.d0*dh*(masRes(i)+4.d0*masRes(i+1)+masRes(i+2))/6.d0
      i=i+2
    end do

    z2Rm=z/Rm
    res=res*8.d0*M*(l2d/50.d0)

    itt=itt+1
    !==after the loop i=2*nn+2, therefore T_new(i) would be out of bounds==!
    hh=find_base(B12,T_new(2*nn+1),alpha,M,R,H,L39,S_d) !; write(*,*)"hh=",hh
  end do
return
end subroutine maxL_n


!============================================================================================
!  Finds the base height h/R where the gas+radiation pressure becomes comparable
!  to the magnetic pressure. The root is found by bisection in the interval [1d-17, 0.1].
!============================================================================================
real*8 function find_base(B12,T,alpha,M,R,H,L39,S_d)
implicit none
real*8,intent(in)::B12,T,alpha,M,R,H,L39,S_d
real*8::P_g_rad
real*8::hh1,hh2,hh
real*8::vel
real*8::columnVel  !==function==!

  hh1=1.d-17
  hh2=0.1d0
  do while(abs(hh1-hh2)/hh2.gt.0.05d0)
    hh=(hh1+hh2)/2.d0
    vel=columnVel(hh*R,H*R,M,R,alpha)
    P_g_rad=7.5d31*L39/S_d/M/vel*1.6d-9*T/0.6d0&
            +5.2d13/3.d0*T**4
    !write(*,*)hh1,hh2,P_g_rad !/B12/B12/1.d24*4.d0*3.141592d0
    if(P_g_rad/B12/B12/1.d24*4.d0*3.1416d0.ge.1.d0)then
      hh1=hh
    else
      hh2=hh
    end if
  end do
  find_base=hh

return
end function find_base
!============================================================================================


!===============================================================================================
!  Diagnostic routine: estimates ram and radiation pressure terms in the accretion channel
!  as a function of magnetic latitude lambda, for given surface field b0 and luminosity L39.
!  Prints: lambda, beta, density, ram pressure, radiation term, rad/ram ratio, local Ecyc.
!===============================================================================================
subroutine channel_pressure(lambda,b0,L39)
implicit none

real*8,intent(in)::b0,L39,lambda

real*8::Rm,S_d,B12,pi
real*8::zA,zB,zC,z,d,length,alpha_disc,M
real*8::ro,M19,ram,rad,beta,b
real*8::cos_lam,cos2_lam,geom_fac,rad_proj,chan_press

  pi=3.14159265359d0

  M=1.4d0
  alpha_disc=0.1d0
  length=0.5d0

  B12=b0*44.13d0
  M19=0.75d0*L39/M
  cos_lam=cos(lambda)
  cos2_lam=cos_lam*cos_lam

  Rm=7.d7*0.5d0*(M**(1.d0/7.d0))*(B12**(4.d0/7.d0))*(L39**(-2.d0/7.d0))
  geom_fac=sqrt((4.d0-3.d0*cos2_lam)/(4.d0-3.d0*(1.d6/Rm)))
  b=b0*((1.d6/Rm)**3)/(cos_lam**6)*geom_fac
  zC=Rm*0.056d0*(L39**(3.d0/20.d0))/(M**(21.d0/40.d0))*((Rm/1.d8)**(1.d0/8.d0))*(alpha_disc**(-0.1d0))     !==C-zone==!
  zB=Rm*0.064d0*(L39**(1.d0/5.d0))/(M**(11.d0/20.d0))*((Rm/1.d8)**(1.d0/20.d0))*(alpha_disc**(-0.1d0))     !==B-zone==!
  zA=1.d7*L39/M                                                                                         !==A-zone==!
  z=max(zA,zB,zC)
  d=1.d6*abs(asin(sqrt(1.d6/Rm))-asin(sqrt(1.d6/(Rm-z))))
  S_d=2.d0*pi*length*1.d6*sqrt(1.d6/(Rm-z))*d

  beta=sqrt(3.d5*M/(Rm*cos2_lam))

  ro=1.d19*M19/2.d0/S_d*((1.d6/Rm)**3)/(cos_lam**5)*geom_fac&
     /beta/3.d10

  ram=ro*(beta**2)*9.d20

  rad=(1.d39*L39-2.d38)/((Rm*cos2_lam)**2)/3.d10/4.d0/pi

  rad_proj=rad*cos(atan(1.d0/2.d0/tan(lambda)))

  chan_press=ram-rad_proj

  !write(*,*)lambda,beta,ro,ram,rad,rad_proj/ram,b*44.13d0*11.6d0
return
end subroutine channel_pressure
!===============================================================================================


!=======================================================================================
!  Finds the maximum value of chan_pressure(lambda,b0,L39)
!  over lambda in the interval [0,1.54] with step 0.01.
!=======================================================================================
real*8 function chan_pres_max(b0,L39)
implicit none
real*8,intent(in)::b0,L39
real*8::chan_pressure  !==function==!
real*8::lambda,res
real*8::lambda_min,lambda_max,dlambda
  lambda_min=0.d0
  lambda_max=1.54d0
  dlambda=0.01d0
  chan_pres_max=-1.d100
  lambda=lambda_min
  do while(lambda.le.lambda_max)
    res=chan_pressure(lambda,b0,L39)
    if(res.gt.chan_pres_max)then
      chan_pres_max=res
    end if
    lambda=lambda+dlambda
  end do
return
end function chan_pres_max
!=======================================================================================


!=======================================================================================
!  Returns the ratio of the projected radiation term to the ram pressure
!  in the accretion channel at magnetic latitude lambda.
!=======================================================================================
real*8 function chan_pressure(lambda,b0,L39)
implicit none
real*8,intent(in)::b0,L39,lambda
real*8::Rm,S_d,B12,pi
real*8::zA,zB,zC,z,d,length,alpha_disc,M
real*8::ro,M19,ram,rad,beta,b
real*8::cos_lam,cos2_lam,tan_lam,geom_fac,proj_fac

  pi=3.14159265359d0

  M=1.4d0
  alpha_disc=0.1d0
  length=0.5d0

  B12=b0*44.13d0
  M19=0.75d0*L39/M

  cos_lam=cos(lambda)
  cos2_lam=cos_lam*cos_lam
  tan_lam=tan(lambda)

  Rm=7.d7*0.5d0*(M**(1.d0/7.d0))*(B12**(4.d0/7.d0))*(L39**(-2.d0/7.d0)) !; write(*,*)"*** ",L39,Rm

  geom_fac=sqrt((4.d0-3.d0*cos2_lam)/(4.d0-3.d0*(1.d6/Rm)))

  b=b0*((1.d6/Rm)**3)/(cos_lam**6)*geom_fac

  zC=Rm*0.056d0*(L39**(3.d0/20.d0))/(M**(21.d0/40.d0))*((Rm/1.d8)**(1.d0/8.d0))*(alpha_disc**(-0.1d0))   !==C-zone==!
  zB=Rm*0.064d0*(L39**(1.d0/5.d0))/(M**(11.d0/20.d0))*((Rm/1.d8)**(1.d0/20.d0))*(alpha_disc**(-0.1d0))   !==B-zone==!
  zA=1.d7*L39/M                                                                                       !==A-zone==!

  z=max(zA,zB,zC)

  d=1.d6*abs(asin(sqrt(1.d6/Rm))-asin(sqrt(1.d6/(Rm-z))))
  S_d=2.d0*pi*length*1.d6*sqrt(1.d6/(Rm-z))*d

  beta=sqrt(3.d5*M/(Rm*cos2_lam))

  ro=1.d19*M19/2.d0/S_d*((1.d6/Rm)**3)/(cos_lam**5)*geom_fac&
     /beta/3.d10

  ram=ro*(beta**2)*9.d20

  rad=(1.d39*L39-2.d38)/((Rm*cos2_lam)**2)/3.d10/4.d0/pi

  proj_fac=2.d0*tan_lam/sqrt(1.d0+4.d0*tan_lam*tan_lam)

  chan_pressure=rad*proj_fac/ram

return
end function chan_pressure
!=======================================================================================





!=============================================================================================
!  Auxiliary analytical integrals and numerical transverse integral used
!  in the accretion-column structure calculations.
!=============================================================================================
real*8 function intF1(x)
implicit none
real*8,intent(in)::x
  intF1=log(x/(x+1.d0))+(12.d0*x**3+42.d0*x**2+52.d0*x+25.d0)/12.d0/(x+1.d0)**4
  !intF1=log((sqrt(1.d0+x)-1.d0)/(sqrt(1.d0+x)+1.d0))+2.d0*(105.d0*(x+1.d0)**3+35.d0*(x+1.d0)**2+21.d0*(x+1.d0)+15.d0)/&
  !                                                   105.d0/((x+1.d0)**3.5d0)       !==Valery correction==!
return
end function intF1
!=================


real*8 function intF2(x)
implicit none
real*8,intent(in)::x
  intF2=5.d0*log((x+1.d0)/x)-(60.d0*x**4+210.d0*x**3+260.d0*x**2+125.d0*x+12.d0)&
                             /12.d0/x/(x+1.d0)**4
return
end function intF2
!================

real*8 function intF3(x)
implicit none
real*8,intent(in)::x
  intF3=15.d0*log(x/(x+1.d0))-(60.d0*x**5+210.d0*x**4+260.d0*x**3+125.d0*x**2+12.d0*x-2.d0)&
                             /4.d0/x/x/(x+1.d0)**4
return
end function intF3
!================

real*8 function intF5(x)
implicit none
real*8,intent(in)::x
  !intF5=70.d0*log(x/(x+1.d0))+(840.d0*x**7+2940.d0*x**6+3640.d0*x**5+1750.d0*x**4+168.d0*x**3-28.d0*x**2+8.d0*x-3.d0)/&
  !      (x+1.d0)**4/12.d0/x**4
  intF5=70.d0*log(x/(x+1.d0))+(840.d0*x**7+2940.d0*x**6+3640.d0*x**5+1750.d0*x**4+168.d0*x**3-28.d0*x**2+8.d0*x-3.d0)/(x+1.d0)**4&
        /12.d0/x**4
return
end function intF5
!=================

real*8 function I1(h,H0,alpha)
implicit none
real*8,intent(in)::h,H0,alpha
real*8::intF1,intF2,intF3,intF5  !==functions==!
integer::a
  a=int(alpha)
  select case(a)
    case(1)
      I1=H0*(intF1(H0)-intF1(h))
    case(2)
      I1=(H0**2)*(intF2(H0)-intF2(h))
    case(3)
      !==old code used intF2 here; check whether intF3 was intended==!
      I1=(H0**3)*(intF2(H0)-intF2(h))
      !I1=(H0**3)*(intF3(H0)-intF3(h))
    case(5)
      I1=(H0**5)*(intF5(H0)-intF5(h))
    case default
      I1=H0*(intF1(H0)-intF1(h))
  end select
return
end function I1
!=========================================================================================================================================

real*8 function I2(b,x,T,n)
implicit none
real*8,intent(in)::b,x
integer,intent(in)::n
real*8,intent(in)::T(n,2)
real*8::sigma_average_ross_ort,find_H  !==functions==!
real*8::y,dy,T_
integer::nn,i
real*8::masI2(13)

  nn=12
  if(x.ge.1.d0)then
    I2=0.d0
  else
    dy=(1.d0-x)/nn
    do i=1,nn+1
      y=x+(i-1)*dy
      T_=find_H(y,T,n)
      !write(*,*)!"t ",T_
      masI2(i)=sigma_average_ross_ort(b,T_)*y
    end do
    I2=0.d0
    do i=1,nn-1,2
      I2=I2+2.d0*dy*(masI2(i)+4.d0*masI2(i+1)+masI2(i+2))/6.d0
    end do
  end if
return
end function I2
!=====================================================================================

!=============================================================================
! The function gives velosity profile in the accretion column
! h - given height [km]; Hcol - accretion column height [km]
! Hmax - column height
! M - neutron star mass measured in solar masses, R - [NS radius}/1[km}
!=============================================================================
real*8 function columnVel(h,Hcol,M,R,alpha)
implicit none
real*8,intent(in)::h,Hcol,M,R,alpha
real*8::Rsh,v_max,v_min
  !alpha=5.d0
  Rsh=2.95d0*M
  v_max=sqrt(Rsh/(R+Hcol))/7.d0
  v_min=0.d0 !1.d-7    !==??
  if(h.le.Hcol)then
    columnVel=v_min+(v_max-v_min)*((h/Hcol)**alpha)
    columnVel=columnVel*sqrt((R+Hcol)/(R+h))     !==Valery correction==!
  else
    columnVel=sqrt(Rsh/(R+Hcol))   !==Just free-fall velocity
  end if
return
end function columnVel
!=============================================================================


!===================================================================================
!  Rosseland mean cross-section for photons propagating along the magnetic field.
!  b - magnetic field strength in units of critical B-field strength
!  T - temperature [keV]
!  Output is in units of the Thomson cross-section.
!===================================================================================
real*8 function sigma_average_ross_par(b,T)
use mod_sig_ave
implicit none
interface
  real*8 function f_sig_ave_ross(x)
    real*8,intent(in)::x
  end function f_sig_ave_ross

  real*8 function f_pl_ave(x)
    real*8,intent(in)::x
  end function f_pl_ave
end interface

real*8,intent(in)::b,T
integer::l_pol_1,l_pol_2
real*8::sigma_ave_1,sigma_ave_2,lim1,lim2,eps
real*8::int_simpson_new  !==function==!
  eps=0.03d0
  !==polarization states; currently pure X-mode if l_pol_1=l_pol_2=1==!
  l_pol_1=1
  l_pol_2=1

  !==we use these variables in module mod_sig_ave==!
  bb=b
  TT=T
  lim1=TT/511.d0/100.d0
  lim2=TT/511.d0*40.d0

  if(l_pol_1.eq.l_pol_2)then
    l_pol=l_pol_1
    sigma_ave_1=1.d0/int_simpson_new(f_sig_ave_ross,lim1,lim2,eps)*&
                int_simpson_new(f_pl_ave,lim1,lim2,eps)
    sigma_ave_2=sigma_ave_1
  else
    l_pol=l_pol_1
    sigma_ave_1=1.d0/int_simpson_new(f_sig_ave_ross,lim1,lim2,eps)*&
                int_simpson_new(f_pl_ave,lim1,lim2,eps)
    l_pol=l_pol_2
    sigma_ave_2=1.d0/int_simpson_new(f_sig_ave_ross,lim1,lim2,eps)*&
                int_simpson_new(f_pl_ave,lim1,lim2,eps)
  end if

  sigma_average_ross_par=1.d0/(0.5d0/sigma_ave_1+0.5d0/sigma_ave_2)
  !write(*,*)sigma_average_ross_par
return
end function sigma_average_ross_par
!================================================================================


!================================================================================
!  Rosseland mean cross-section for photons propagating across the magnetic field.
!  b - magnetic field strength in units of critical B-field strength
!  T - temperature [keV]
!  Output is in units of the Thomson cross-section.
!===============================================================================
real*8 function sigma_average_ross_ort(b,T)
use mod_sig_ave
implicit none
interface
  real*8 function f_sig_ave_ross_ort(x)
    real*8,intent(in)::x
  end function f_sig_ave_ross_ort

  real*8 function f_pl_ave(x)
    real*8,intent(in)::x
  end function f_pl_ave
end interface

real*8,intent(in)::b,T
integer::l_pol_1,l_pol_2
real*8::sigma_ave_1,sigma_ave_2,lim1,lim2,eps
real*8::int_simpson_new  !==function==!

  eps=0.03d0

  !==polarization states; currently pure X-mode if l_pol_1=l_pol_2=1==!
  l_pol_1=1
  l_pol_2=1

  !==we use these variables in module mod_sig_ave==!
  bb=b
  TT=T
  lim1=TT/511.d0/100.d0
  lim2=TT/511.d0*40.d0

  if(l_pol_1.eq.l_pol_2)then
    l_pol=l_pol_1
    sigma_ave_1=1.d0/int_simpson_new(f_sig_ave_ross_ort,lim1,lim2,eps)*&
                int_simpson_new(f_pl_ave,lim1,lim2,eps)
    sigma_ave_2=sigma_ave_1
  else
    l_pol=l_pol_1
    sigma_ave_1=1.d0/int_simpson_new(f_sig_ave_ross_ort,lim1,lim2,eps)*&
                int_simpson_new(f_pl_ave,lim1,lim2,eps)
    l_pol=l_pol_2
    sigma_ave_2=1.d0/int_simpson_new(f_sig_ave_ross_ort,lim1,lim2,eps)*&
                int_simpson_new(f_pl_ave,lim1,lim2,eps)

  end if
  sigma_average_ross_ort=1.d0/(0.5d0/sigma_ave_1+0.5d0/sigma_ave_2)
return
end function sigma_average_ross_ort
!================================================================================


!================================================================================================================
!  Rosseland integrand for photons propagating across the magnetic field.
!  x - photon energy in units of electron rest energy
!  TT, bb, l_pol and the cross-section tables are taken from mod_sig_ave.
!================================================================================================================
real*8 function f_sig_ave_ross_ort(x)
use mod_sig_ave
implicit none
real*8,intent(in)::x
real*8::pi
real*8::sig_ang_ave_ross_ort,dplank_dT  !==functions==!
  pi=3.1415926535d0
  f_sig_ave_ross_ort=6.d0/pi*dplank_dT(x*511.d3*1.602d-12,TT*1.d3*11600.d0)&
                     /sig_ang_ave_ross_ort(x,bb,sig1,sig2,10,l_pol)
return
end function f_sig_ave_ross_ort
!================================================================================================================


!==angle-average-Rosseland CS=================================================
! k - photon energy
! b - magnetic field strength
! sig1 and sig2 - arrays with the scattering cross sections for X- and O-modes
! n - number of angles for averaging
! l - polarization state
!=============================================================================
real*8 function sig_ang_ave_ross_ort(k,b,sig1,sig2,n,l)
implicit none
integer,intent(in)::n,l
real*8,intent(in)::k,b,sig1,sig2
dimension sig1(n,300,2), sig2(n,300,2)
real*8::mas
dimension mas(300,2)
real*8::pi,dthetta,thetta,dfi,fi,thettaB,s1,s2,s
integer::i
real*8::find_H,cross_sect_app,int_magn_knt  !==functions==!
  pi=3.1415926535d0

  dthetta=pi/2/n
  thetta=dthetta/2
  i=1
  dfi=pi/2/n
  sig_ang_ave_ross_ort=0.d0
  do while(thetta.le.pi/2)
    fi=dfi/2
    thettaB=acos(sin(thetta)*cos(fi))
    i=thettaB/dthetta
    i=max(1,i)
    do while(fi.le.pi/2.d0)
      if(k.gt.(b/10.d0))then
        if(k.lt.(3.2d0*b))then
          mas(1:300,1) = sig1(i,1:300,1); mas(1:300,2) = sig1(i,1:300,2)
          s1=int_magn_knt(k,thetta,b,1)
          !s1=find_H(k,mas,300)
          mas(1:300,1) = sig2(i,1:300,1); mas(1:300,2) = sig2(i,1:300,2)
          s2=int_magn_knt(k,thetta,b,2)
          !s2=find_H(k,mas,300)
        else
          s1 = int_magn_knt(k,thetta,b,1); s2 = int_magn_knt(k,thetta,b,2)
        end if
      else
        s1 = int_magn_knt(k,thetta,b,1)
        s2 = int_magn_knt(k,thetta,b,2)
        !s1=1.d0; s2=1.d0
      end if
      if(l.eq.1)then
        s = s1
      else
        s = s2
      end if
      !s=1.d0
      sig_ang_ave_ross_ort = sig_ang_ave_ross_ort+dfi*dthetta*sin(thetta)*cos(thetta)*cos(thetta)/s
      fi = fi+dfi
    end do
    thetta=thetta+dthetta
  end do
  sig_ang_ave_ross_ort=1.d0/sig_ang_ave_ross_ort
return
end function sig_ang_ave_ross_ort
!=======================================================================================


!=======================================================================
!  Derivative of the Planck function with respect to temperature.
!  [E]=[erg], [T]=[K]
!  [dplank_dT]=dB_nu/dT in cgs units
!=======================================================================
real*8 function dplank_dT(E,T)
implicit none
real*8,intent(in)::E,T
real*8::c,h,k
real*8::u,exp_u
  c=2.99792d10
  h=6.62606896d-27
  k=1.3807d-16
  u=E/(k*T)
  if(u.gt.100.d0)then
    dplank_dT=2.d0*(E**3)/((h**2)*(c**2))*exp(-u)*E/k/(T**2)
  else if(u.lt.0.01d0)then
    dplank_dT=2.d0*k*(E**2)/((h**2)*(c**2))
  else
    exp_u=exp(u)
    dplank_dT=2.d0*(E**3)/((h**2)*(c**2))&
              *E*exp_u/k/(T**2)/(exp_u-1.d0)**2
  end if
return
end function dplank_dT
!=======================================================================


!=======================================================================
!  Planck weight for Rosseland averaging.
!  x - photon energy in units of electron rest energy.
!  TT is temperature [keV] from mod_sig_ave.
!=======================================================================
real*8 function f_pl_ave(x)
use mod_sig_ave
implicit none
real*8,intent(in)::x
real*8::dplank_dT  !==function==!
  !f_pl_ave=plank_qft_norn(x,TT)
  f_pl_ave=dplank_dT(x*511.d3*1.602d-12,TT*1.d3*11600.d0)
return
end function f_pl_ave
!=======================================================================


!================================================================================================================
!  Rosseland integrand for photons propagating along the magnetic field.
!  x - photon energy in units of electron rest energy.
!  TT, bb, l_pol and the cross-section tables are taken from mod_sig_ave.
!================================================================================================================
real*8 function f_sig_ave_ross(x)
use mod_sig_ave
implicit none
real*8,intent(in)::x
real*8::sig_ang_ave_ross,dplank_dT  !==functions==!
  f_sig_ave_ross=3.d0*dplank_dT(x*511.d3*1.602d-12,TT*1.d3*11600.d0)&
                 /sig_ang_ave_ross(x,bb,sig1,sig2,10,l_pol)
  !write(*,*)x,sig_ang_ave_ross(x,bb,sig1,sig2,10,l_pol); read(*,*)
  !f_sig_ave_ross=dplank_dT(x*511.d3*1.602d-12,TT*1.d3*11600.d0)/((x/bb)**2)
return
end function f_sig_ave_ross
!================================================================================================================


!====================================================================================
!  Approximate magnetic cross-section without resonances.
!  k      - photon energy in units of electron rest energy
!  thetta_ - angle between photon momentum and magnetic field
!  b      - magnetic field strength in units of critical field
!  l      - polarization state
!====================================================================================
real*8 function int_magn_knt(k,thetta_,b,l)
implicit none
real*8,intent(in)::k,thetta_,b
integer,intent(in)::l
real*8::k_res
real*8::sin_t
real*8::int_cs_knt_app   !==function==!
  sin_t=sin(thetta_)
  if(thetta_.eq.0.d0)then
    k_res=b
  else
    !k_res=b !(sqrt(1.d0+2.d0*b*(sin_t**2))-1.d0)/(sin_t**2)  !==check this approximation==!
    k_res=(sqrt(1.d0+2.d0*b*(sin_t**2))-1.d0)/(sin_t**2)  !==check this approximation==!
  end if
  if(l.eq.1)then
    int_magn_knt=min(int_cs_knt_app(k),(k/k_res)**2)
  else
    int_magn_knt=min(int_cs_knt_app(k),(k/k_res)**2+sin_t**2)
  end if
return
end function int_magn_knt
!====================================================================================

!==========================================================================================
!==========================================================================================
real*8 function int_cs_knt_app(k)
implicit none
real*8,intent(in)::k
real*8::knt
dimension knt(72,2)
real*8::find_H  !==function==!
  !====================================================
   knt(           1 ,1)=   1.0000000000000000E-004 ;  knt(           1 ,2)=   1.0000973202057453
   knt(           2 ,1)=   1.2000000000000000E-004 ;  knt(           2 ,2)=   1.0000573311085399
   knt(           3 ,1)=   1.4400000000000000E-004 ;  knt(           3 ,2)=   1.0000093496805458
   knt(           4 ,1)=   1.7280000000000000E-004 ;  knt(           4 ,2)=  0.99995177986547124
   knt(           5 ,1)=   2.0735999999999999E-004 ;  knt(           5 ,2)=  0.99988270746308550
   knt(           6 ,1)=   2.4883199999999999E-004 ;  knt(           6 ,2)=  0.99979983695689978
   knt(           7 ,1)=   2.9859839999999999E-004 ;  knt(           7 ,2)=  0.99970041591930492
   knt(           8 ,1)=   3.5831808000000000E-004 ;  knt(           8 ,2)=  0.99958114460833258
   knt(           9 ,1)=   4.2998169600000000E-004 ;  knt(           9 ,2)=  0.99943806787219491
   knt(          10 ,1)=   5.1597803519999998E-004 ;  knt(          10 ,2)=  0.99926644607903836
   knt(          11 ,1)=   6.1917364223999999E-004 ;  knt(          11 ,2)=  0.99906060107692773
   knt(          12 ,1)=   7.4300837068800002E-004 ;  knt(          12 ,2)=  0.99881373260860673
   knt(          13 ,1)=   8.9161004482559997E-004 ;  knt(          13 ,2)=  0.99851769982289973
   knt(          14 ,1)=   1.0699320537907199E-003 ;  knt(          14 ,2)=  0.99816276162777928
   knt(          15 ,1)=   1.2839184645488638E-003 ;  knt(          15 ,2)=  0.99773726885356340
   knt(          16 ,1)=   1.5407021574586365E-003 ;  knt(          16 ,2)=  0.99722730009714178
   knt(          17 ,1)=   1.8488425889503638E-003 ;  knt(          17 ,2)=  0.99661623231366347
   knt(          18 ,1)=   2.2186111067404365E-003 ;  knt(          18 ,2)=  0.99588423631526024
   knt(          19 ,1)=   2.6623333280885236E-003 ;  knt(          19 ,2)=  0.99500768673016127
   knt(          20 ,1)=   3.1947999937062283E-003 ;  knt(          20 ,2)=  0.99395847582230035
   knt(          21 ,1)=   3.8337599924474740E-003 ;  knt(          21 ,2)=  0.99270322109376585
   knt(          22 ,1)=   4.6005119909369686E-003 ;  knt(          22 ,2)=  0.99120235821964298
   knt(          23 ,1)=   5.5206143891243621E-003 ;  knt(          23 ,2)=  0.98940911434654710
   knt(          24 ,1)=   6.6247372669492348E-003 ;  knt(          24 ,2)=  0.98726836265574580
   knt(          25 ,1)=   7.9496847203390821E-003 ;  knt(          25 ,2)=  0.98471536859449449
   knt(          26 ,1)=   9.5396216644068974E-003 ;  knt(          26 ,2)=  0.98167445266317044
   knt(          27 ,1)=   1.1447545997288276E-002 ;  knt(          27 ,2)=  0.97805761565593607
   knt(          28 ,1)=   1.3737055196745932E-002 ;  knt(          28 ,2)=  0.97376320162824004
   knt(          29 ,1)=   1.6484466236095119E-002 ;  knt(          29 ,2)=  0.96867471362699786
   knt(          30 ,1)=   1.9781359483314141E-002 ;  knt(          30 ,2)=  0.96265994836396851
   knt(          31 ,1)=   2.3737631379976969E-002 ;  knt(          31 ,2)=  0.95557067909315563
   knt(          32 ,1)=   2.8485157655972360E-002 ;  knt(          32 ,2)=  0.94724318798818741
   knt(          33 ,1)=   3.4182189187166832E-002 ;  knt(          33 ,2)=  0.93750002374799690
   knt(          34 ,1)=   4.1018627024600199E-002 ;  knt(          34 ,2)=  0.92615342294705283
   knt(          35 ,1)=   4.9222352429520236E-002 ;  knt(          35 ,2)=  0.91301086150540756
   knt(          36 ,1)=   5.9066822915424283E-002 ;  knt(          36 ,2)=  0.89788316158320425
   knt(          37 ,1)=   7.0880187498509134E-002 ;  knt(          37 ,2)=  0.88059542506014443
   knt(          38 ,1)=   8.5056224998210958E-002 ;  knt(          38 ,2)=  0.86100075098434048
   knt(          39 ,1)=  0.10206746999785314      ;  knt(          39 ,2)=  0.83899618864718462
   knt(          40 ,1)=  0.12248096399742377      ;  knt(          40 ,2)=  0.81453968798130294
   knt(          41 ,1)=  0.14697715679690851      ;  knt(          41 ,2)=  0.78766601848973361
   knt(          42 ,1)=  0.17637258815629020      ;  knt(          42 ,2)=  0.75849892308673672
   knt(          43 ,1)=  0.21164710578754822      ;  knt(          43 ,2)=  0.72725644002748724
   knt(          44 ,1)=  0.25397652694505785      ;  knt(          44 ,2)=  0.69424668661770439
   knt(          45 ,1)=  0.30477183233406940      ;  knt(          45 ,2)=  0.65985267323616448
   knt(          46 ,1)=  0.36572619880088325      ;  knt(          46 ,2)=  0.62450684879351470
   knt(          47 ,1)=  0.43887143856105987      ;  knt(          47 ,2)=  0.58865861143133635
   knt(          48 ,1)=  0.52664572627327177      ;  knt(          48 ,2)=  0.55274014351938505
   knt(          49 ,1)=  0.63197487152792609      ;  knt(          49 ,2)=  0.51713677678331205
   knt(          50 ,1)=  0.75836984583351130      ;  knt(          50 ,2)=  0.48216715467534654
   knt(          51 ,1)=  0.91004381500021347      ;  knt(          51 ,2)=  0.44807592256416745
   knt(          52 ,1)=   1.0920525780002561      ;  knt(          52 ,2)=  0.41503841520732532
   knt(          53 ,1)=   1.3104630936003072      ;  knt(          53 ,2)=  0.38317399328589180
   knt(          54 ,1)=   1.5725557123203686      ;  knt(          54 ,2)=  0.35256321530585011
   knt(          55 ,1)=   1.8870668547844422      ;  knt(          55 ,2)=  0.32326418713517624
   knt(          56 ,1)=   2.2644802257413308      ;  knt(          56 ,2)=  0.29532482555068029
   knt(          57 ,1)=   2.7173762708895968      ;  knt(          57 ,2)=  0.26878966992995967
   knt(          58 ,1)=   3.2608515250675159      ;  knt(          58 ,2)=  0.24370157390552352
   knt(          59 ,1)=   3.9130218300810187      ;  knt(          59 ,2)=  0.22009968747491873
   knt(          60 ,1)=   4.6956261960972219      ;  knt(          60 ,2)=  0.19801551716126464
   knt(          61 ,1)=   5.6347514353166659      ;  knt(          61 ,2)=  0.17746868171513352
   knt(          62 ,1)=   6.7617017223799989      ;  knt(          62 ,2)=  0.15846351150818441
   knt(          63 ,1)=   8.1140420668559976      ;  knt(          63 ,2)=  0.14098709761782055
   knt(          64 ,1)=   9.7368504802271971      ;  knt(          64 ,2)=  0.12500893146362141
   knt(          65 ,1)=   11.684220576272637      ;  knt(          65 ,2)=  0.11048195481673427
   knt(          66 ,1)=   14.021064691527164      ;  knt(          66 ,2)=   9.7344668409412582E-002
   knt(          67 ,1)=   16.825277629832595      ;  knt(          67 ,2)=   8.5523896488402182E-002
   knt(          68 ,1)=   20.190333155799113      ;  knt(          68 ,2)=   7.4937834137906237E-002
   knt(          69 ,1)=   24.228399786958935      ;  knt(          69 ,2)=   6.5499076512244694E-002
   knt(          70 ,1)=   29.074079744350719      ;  knt(          70 ,2)=   5.7117415673771856E-002
   knt(          71 ,1)=   34.888895693220860      ;  knt(          71 ,2)=   4.9702273235762852E-002
   knt(          72 ,1)=   41.866674831865033      ;  knt(          72 ,2)=   4.3164705807852290E-002
  !=======================================================================================================
  if(k.le.1.d-4)then
    int_cs_knt_app=1.d0
  else
    if(k.le.41.d0)then
      int_cs_knt_app=find_H(k,knt,72)
    else
      int_cs_knt_app=0.75d0/2.d0/k*(log(2.d0*k)+0.5d0)
    end if
  end if
return
end function int_cs_knt_app


!==angle-average-Rosseland CS==!
real*8 function sig_ang_ave_ross(k,b,sig1,sig2,n,l)
implicit none
integer,intent(in)::n,l
real*8,intent(in)::k,b,sig1,sig2
dimension sig1(n,300,2), sig2(n,300,2)
real*8::mas,knt
dimension mas(300,2),knt(72,2)
real*8::pi,dthetta,thetta,s1,s2,s
integer::i
real*8::find_H,cross_sect_app,int_magn_knt  !==functions==!
  pi=3.1415926535d0

  !====================================================
   knt(           1 ,1)=   1.0000000000000000E-004 ;  knt(           1 ,2)=   1.0000973202057453
   knt(           2 ,1)=   1.2000000000000000E-004 ;  knt(           2 ,2)=   1.0000573311085399
   knt(           3 ,1)=   1.4400000000000000E-004 ;  knt(           3 ,2)=   1.0000093496805458
   knt(           4 ,1)=   1.7280000000000000E-004 ;  knt(           4 ,2)=  0.99995177986547124
   knt(           5 ,1)=   2.0735999999999999E-004 ;  knt(           5 ,2)=  0.99988270746308550
   knt(           6 ,1)=   2.4883199999999999E-004 ;  knt(           6 ,2)=  0.99979983695689978
   knt(           7 ,1)=   2.9859839999999999E-004 ;  knt(           7 ,2)=  0.99970041591930492
   knt(           8 ,1)=   3.5831808000000000E-004 ;  knt(           8 ,2)=  0.99958114460833258
   knt(           9 ,1)=   4.2998169600000000E-004 ;  knt(           9 ,2)=  0.99943806787219491
   knt(          10 ,1)=   5.1597803519999998E-004 ;  knt(          10 ,2)=  0.99926644607903836
   knt(          11 ,1)=   6.1917364223999999E-004 ;  knt(          11 ,2)=  0.99906060107692773
   knt(          12 ,1)=   7.4300837068800002E-004 ;  knt(          12 ,2)=  0.99881373260860673
   knt(          13 ,1)=   8.9161004482559997E-004 ;  knt(          13 ,2)=  0.99851769982289973
   knt(          14 ,1)=   1.0699320537907199E-003 ;  knt(          14 ,2)=  0.99816276162777928
   knt(          15 ,1)=   1.2839184645488638E-003 ;  knt(          15 ,2)=  0.99773726885356340
   knt(          16 ,1)=   1.5407021574586365E-003 ;  knt(          16 ,2)=  0.99722730009714178
   knt(          17 ,1)=   1.8488425889503638E-003 ;  knt(          17 ,2)=  0.99661623231366347
   knt(          18 ,1)=   2.2186111067404365E-003 ;  knt(          18 ,2)=  0.99588423631526024
   knt(          19 ,1)=   2.6623333280885236E-003 ;  knt(          19 ,2)=  0.99500768673016127
   knt(          20 ,1)=   3.1947999937062283E-003 ;  knt(          20 ,2)=  0.99395847582230035
   knt(          21 ,1)=   3.8337599924474740E-003 ;  knt(          21 ,2)=  0.99270322109376585
   knt(          22 ,1)=   4.6005119909369686E-003 ;  knt(          22 ,2)=  0.99120235821964298
   knt(          23 ,1)=   5.5206143891243621E-003 ;  knt(          23 ,2)=  0.98940911434654710
   knt(          24 ,1)=   6.6247372669492348E-003 ;  knt(          24 ,2)=  0.98726836265574580
   knt(          25 ,1)=   7.9496847203390821E-003 ;  knt(          25 ,2)=  0.98471536859449449
   knt(          26 ,1)=   9.5396216644068974E-003 ;  knt(          26 ,2)=  0.98167445266317044
   knt(          27 ,1)=   1.1447545997288276E-002 ;  knt(          27 ,2)=  0.97805761565593607
   knt(          28 ,1)=   1.3737055196745932E-002 ;  knt(          28 ,2)=  0.97376320162824004
   knt(          29 ,1)=   1.6484466236095119E-002 ;  knt(          29 ,2)=  0.96867471362699786
   knt(          30 ,1)=   1.9781359483314141E-002 ;  knt(          30 ,2)=  0.96265994836396851
   knt(          31 ,1)=   2.3737631379976969E-002 ;  knt(          31 ,2)=  0.95557067909315563
   knt(          32 ,1)=   2.8485157655972360E-002 ;  knt(          32 ,2)=  0.94724318798818741
   knt(          33 ,1)=   3.4182189187166832E-002 ;  knt(          33 ,2)=  0.93750002374799690
   knt(          34 ,1)=   4.1018627024600199E-002 ;  knt(          34 ,2)=  0.92615342294705283
   knt(          35 ,1)=   4.9222352429520236E-002 ;  knt(          35 ,2)=  0.91301086150540756
   knt(          36 ,1)=   5.9066822915424283E-002 ;  knt(          36 ,2)=  0.89788316158320425
   knt(          37 ,1)=   7.0880187498509134E-002 ;  knt(          37 ,2)=  0.88059542506014443
   knt(          38 ,1)=   8.5056224998210958E-002 ;  knt(          38 ,2)=  0.86100075098434048
   knt(          39 ,1)=  0.10206746999785314      ;  knt(          39 ,2)=  0.83899618864718462
   knt(          40 ,1)=  0.12248096399742377      ;  knt(          40 ,2)=  0.81453968798130294
   knt(          41 ,1)=  0.14697715679690851      ;  knt(          41 ,2)=  0.78766601848973361
   knt(          42 ,1)=  0.17637258815629020      ;  knt(          42 ,2)=  0.75849892308673672
   knt(          43 ,1)=  0.21164710578754822      ;  knt(          43 ,2)=  0.72725644002748724
   knt(          44 ,1)=  0.25397652694505785      ;  knt(          44 ,2)=  0.69424668661770439
   knt(          45 ,1)=  0.30477183233406940      ;  knt(          45 ,2)=  0.65985267323616448
   knt(          46 ,1)=  0.36572619880088325      ;  knt(          46 ,2)=  0.62450684879351470
   knt(          47 ,1)=  0.43887143856105987      ;  knt(          47 ,2)=  0.58865861143133635
   knt(          48 ,1)=  0.52664572627327177      ;  knt(          48 ,2)=  0.55274014351938505
   knt(          49 ,1)=  0.63197487152792609      ;  knt(          49 ,2)=  0.51713677678331205
   knt(          50 ,1)=  0.75836984583351130      ;  knt(          50 ,2)=  0.48216715467534654
   knt(          51 ,1)=  0.91004381500021347      ;  knt(          51 ,2)=  0.44807592256416745
   knt(          52 ,1)=   1.0920525780002561      ;  knt(          52 ,2)=  0.41503841520732532
   knt(          53 ,1)=   1.3104630936003072      ;  knt(          53 ,2)=  0.38317399328589180
   knt(          54 ,1)=   1.5725557123203686      ;  knt(          54 ,2)=  0.35256321530585011
   knt(          55 ,1)=   1.8870668547844422      ;  knt(          55 ,2)=  0.32326418713517624
   knt(          56 ,1)=   2.2644802257413308      ;  knt(          56 ,2)=  0.29532482555068029
   knt(          57 ,1)=   2.7173762708895968      ;  knt(          57 ,2)=  0.26878966992995967
   knt(          58 ,1)=   3.2608515250675159      ;  knt(          58 ,2)=  0.24370157390552352
   knt(          59 ,1)=   3.9130218300810187      ;  knt(          59 ,2)=  0.22009968747491873
   knt(          60 ,1)=   4.6956261960972219      ;  knt(          60 ,2)=  0.19801551716126464
   knt(          61 ,1)=   5.6347514353166659      ;  knt(          61 ,2)=  0.17746868171513352
   knt(          62 ,1)=   6.7617017223799989      ;  knt(          62 ,2)=  0.15846351150818441
   knt(          63 ,1)=   8.1140420668559976      ;  knt(          63 ,2)=  0.14098709761782055
   knt(          64 ,1)=   9.7368504802271971      ;  knt(          64 ,2)=  0.12500893146362141
   knt(          65 ,1)=   11.684220576272637      ;  knt(          65 ,2)=  0.11048195481673427
   knt(          66 ,1)=   14.021064691527164      ;  knt(          66 ,2)=   9.7344668409412582E-002
   knt(          67 ,1)=   16.825277629832595      ;  knt(          67 ,2)=   8.5523896488402182E-002
   knt(          68 ,1)=   20.190333155799113      ;  knt(          68 ,2)=   7.4937834137906237E-002
   knt(          69 ,1)=   24.228399786958935      ;  knt(          69 ,2)=   6.5499076512244694E-002
   knt(          70 ,1)=   29.074079744350719      ;  knt(          70 ,2)=   5.7117415673771856E-002
   knt(          71 ,1)=   34.888895693220860      ;  knt(          71 ,2)=   4.9702273235762852E-002
   knt(          72 ,1)=   41.866674831865033      ;  knt(          72 ,2)=   4.3164705807852290E-002
  !=======================================================================================================


  dthetta=pi/2.d0/n; thetta=dthetta/2.d0; i=1
  sig_ang_ave_ross=0.d0
  do while(thetta.le.pi/2.d0)
    if(k.gt.(b/10.d0))then
      if(k.lt.(3.2d0*b))then
        mas(1:300,1)=sig1(i,1:300,1); mas(1:300,2)=sig1(i,1:300,2)
        !s1=int_magn_knt(k,thetta,0.d0,b,1)
        s1=find_H(k,mas,300)
        mas(1:300,1)=sig2(i,1:300,1); mas(1:300,2)=sig2(i,1:300,2)
        !s2=int_magn_knt(k,thetta,0.d0,b,2)
        s2=find_H(k,mas,300)
      else
        if(k.lt.40.d0)then
          s1=find_H(k,knt,72); s2=s1
        else
          s1=0.75d0/2.d0/k*(log(2.d0*k)+0.5d0);
        end if
        !s1=cross_sect_app(k,b); s2=s1   !KNT-tail
        !s1=1.d0;  s2=s1
      end if
    else
      !s1=((k/b)**2)*(1.d0+(cos(thetta))**2); s2=(sin(thetta))**2
      s1=((k/b)**2); s2=(sin(thetta))**2+((cos(thetta))**2)*((k/b)**2)
      !s1=1.d0; s2=1.d0
    end if
    if(l.eq.1)then
      s=s1
    else
      s=s2
    end if
    !s=1.d0
    sig_ang_ave_ross=sig_ang_ave_ross+dthetta*sin(thetta)*cos(thetta)*cos(thetta)/s
    !/((exp(1.d0)-1.d0)/(exp(1.d0/1.1d0/(1.d0+0.4d0*cos(thetta)))-1.d0))
    thetta=thetta+dthetta; i=i+1
    !write(*,*)thetta
  end do
  sig_ang_ave_ross=1.d0/sig_ang_ave_ross
return
end function sig_ang_ave_ross
!=======================================================================================


!===========================================================================
!===========================================================================
subroutine test_AC()
use omp_lib
use iso_fortran_env, only: int64
implicit none
integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: G        = 6.67430d-8
real(dp), parameter :: c_light  = 2.99792458d10
real(dp), parameter :: sigma_T  = 6.6524587321d-25
real(dp), parameter :: m_p      = 1.67262192369d-24
real(dp), parameter :: k_B      = 1.380649d-16
real(dp), parameter :: a_rad    = 7.5657d-15
real(dp), parameter :: M_sun    = 1.98847d33
real(dp), parameter :: erg_keV  = 1.602176634d-9
real(dp), parameter :: pi       = 3.1415926535897932384626433832795d0

real(dp), parameter :: M_NS = 1.4d0 * M_sun
real(dp), parameter :: R_NS = 1.0d6

real(dp), parameter :: R_COL = 1.0d4
real(dp), parameter :: H_COL = 3.0d5

real(dp), parameter :: MDOT_COL = 5.d17

integer, parameter :: NR = 40
integer, parameter :: NZ = 120

integer, parameter :: N_PACKETS = 10000 !30000 !15000
integer, parameter :: N_ITER    = 200

real(dp), parameter :: OPACITY_SCALE = 1.0d0

integer, parameter :: MAX_SCATTERS = 1000000
integer, parameter :: MAX_STEPS    = 1000000

real(dp), parameter :: V_MIN = 1.0d6
real(dp), parameter :: VELOCITY_RELAX = 0.07d0 !0.07d0
real(dp), parameter :: TEMP_AVG_ALPHA = 0.3d0

integer, parameter :: N_STREAMS = 20

real(dp) :: r_edges(0:NR), z_edges(0:NZ)
real(dp) :: r_cent(NR), z_cent(NZ)
real(dp) :: dr, dz

real(dp) :: cell_vol(NR,NZ)
real(dp) :: ring_area(NR)
real(dp) :: side_area_zbin

real(dp) :: R2D(NR,NZ), Z2D(NR,NZ)
real(dp) :: g_grid(NR,NZ)
real(dp) :: v_ff_profile(NZ)
real(dp) :: v_ff_grid(NR,NZ)

real(dp) :: v_grid(NR,NZ), v_new(NR,NZ), v_target(NR,NZ)
real(dp) :: rho_grid(NR,NZ), ne_grid(NR,NZ), alpha_grid(NR,NZ)
real(dp) :: dv_dz_grid(NR,NZ)

real(dp) :: q_grid(NR,NZ), lum_cells(NR,NZ)
real(dp) :: u_accum(NR,NZ), Fr_accum(NR,NZ), Fz_accum(NR,NZ)
real(dp) :: u_rad(NR,NZ), Fr(NR,NZ), Fz(NR,NZ)
real(dp) :: u_rad_avg(NR,NZ)
real(dp) :: u_rad_used(NR,NZ)
logical :: have_u_rad_avg
real(dp) :: T_K(NR,NZ), T_keV(NR,NZ)
real(dp) :: P_gas(NR,NZ), P_rad(NR,NZ)
real(dp) :: dPgas_dz(NR,NZ), dPrad_dz(NR,NZ)
real(dp) :: dv_dz_suggested(NR,NZ)

real(dp), allocatable :: u_thr(:,:,:), Fr_thr(:,:,:), Fz_thr(:,:,:)
real(dp), allocatable :: visit_thr(:,:,:), side_thr(:,:)

real(dp) :: cell_cdf(NR*NZ)
real(dp) :: base_ring_cdf(NR)
real(dp) :: base_lum_rings(NR)

real(dp) :: L_volume, L_base, L_total
real(dp) :: packet_luminosity, p_base
real(dp) :: tau_radial_mean, tau_vertical_mean

real(dp) :: area_col
real(dp) :: rho_base(NR), v_base(NR), mdot_ring(NR)

integer :: i, j, k, iteration
integer :: i_vprof

integer(int64) :: seeds(N_STREAMS)

integer :: scatter_hist(0:MAX_SCATTERS+1)
integer, allocatable :: scatter_thr(:,:)

real(dp) :: dv_dz_residual(NR,NZ)
real(dp) :: iter_metric
real(dp), parameter :: RES_EPS = 1.0d-30
real(dp) :: dv_raw, dv_lim

call execute_command_line("mkdir -p ./res")
call execute_command_line("mkdir -p ./res/AC_gif")

call omp_set_num_threads(N_STREAMS)

scatter_hist = 0

do i = 0, NR
   r_edges(i) = R_COL * real(i,dp) / real(NR,dp)
end do

do j = 0, NZ
   z_edges(j) = H_COL * real(j,dp) / real(NZ,dp)
end do

do i = 1, NR
   r_cent(i) = 0.5d0 * (r_edges(i-1) + r_edges(i))
end do

do j = 1, NZ
   z_cent(j) = 0.5d0 * (z_edges(j-1) + z_edges(j))
end do

dr = r_edges(1) - r_edges(0)
dz = z_edges(1) - z_edges(0)

do i = 1, NR
   ring_area(i) = pi * (r_edges(i)**2 - r_edges(i-1)**2)
end do

do i = 1, NR
   do j = 1, NZ
      cell_vol(i,j) = ring_area(i) * dz
      R2D(i,j) = r_cent(i)
      Z2D(i,j) = z_cent(j)
   end do
end do

side_area_zbin = 2.0d0 * pi * R_COL * dz

i_vprof = 1
do i = 1, NR
   if (abs(r_cent(i) - R_COL/3.0d0) < abs(r_cent(i_vprof) - R_COL/3.0d0)) then
      i_vprof = i
   end if
end do

do j = 1, NZ
   v_ff_profile(j) = v_freefall(z_cent(j))
end do

do i = 1, NR
   do j = 1, NZ
      g_grid(i,j) = gravity(z_cent(j))
      v_ff_grid(i,j) = v_ff_profile(j)
      v_grid(i,j) = v_ff_grid(i,j)
   end do
end do

do i = 1, N_STREAMS
   seeds(i) = int(12345 + 104729*i, int64)
end do

allocate(u_thr(NR,NZ,N_STREAMS))
allocate(Fr_thr(NR,NZ,N_STREAMS))
allocate(Fz_thr(NR,NZ,N_STREAMS))
allocate(visit_thr(NR,NZ,N_STREAMS))
allocate(side_thr(NZ,N_STREAMS))
allocate(scatter_thr(0:MAX_SCATTERS+1,N_STREAMS))

dv_dz_residual = 0.0d0
call write_velocity_profile(0)
have_u_rad_avg = .false.
u_rad_avg = 0.0d0
u_rad_used = 0.0d0

do iteration = 1, N_ITER

   write(*,*)
   write(*,*) "================================================"
   write(*,'(A,I4,A,I4)') "ITERATION ", iteration, " / ", N_ITER
   write(*,*) "================================================"

   area_col = pi * R_COL**2

   do i = 1, NR
      do j = 1, NZ
         rho_grid(i,j)   = MDOT_COL / (area_col * v_grid(i,j))
         ne_grid(i,j)    = rho_grid(i,j) / m_p
         alpha_grid(i,j) = OPACITY_SCALE * ne_grid(i,j) * sigma_T
      end do
   end do

   tau_radial_mean   = sum(alpha_grid) / real(NR*NZ,dp) * R_COL
   tau_vertical_mean = sum(alpha_grid) / real(NR*NZ,dp) * H_COL

   write(*,'(A,ES12.4)') "Mean radial optical depth   = ", tau_radial_mean
   write(*,'(A,ES12.4)') "Mean vertical optical depth = ", tau_vertical_mean

   call gradient_z(v_grid, dv_dz_grid, dz)

   do i = 1, NR
      do j = 1, NZ
         q_grid(i,j) = rho_grid(i,j) * v_grid(i,j) * &
                       (g_grid(i,j) + v_grid(i,j) * dv_dz_grid(i,j))
         if (q_grid(i,j) < 0.0d0) q_grid(i,j) = 0.0d0
         lum_cells(i,j) = q_grid(i,j) * cell_vol(i,j)
      end do
   end do

   L_volume = sum(lum_cells)

   do i = 1, NR
      v_base(i) = v_grid(i,1)
      rho_base(i) = rho_grid(i,1)
      mdot_ring(i) = rho_base(i) * v_base(i) * ring_area(i)
      base_lum_rings(i) = 0.5d0 * mdot_ring(i) * v_base(i)**2
   end do

   L_base = sum(base_lum_rings)
   L_total = L_volume + L_base

   write(*,'(A,ES12.4)') "L_volume = ", L_volume
   write(*,'(A,ES12.4)') "L_base   = ", L_base
   write(*,'(A,ES12.4)') "L_total  = ", L_total

   if (L_total <= 0.0d0) stop "ERROR: total luminosity is non-positive."

   call build_cdf_2d(lum_cells, NR, NZ, cell_cdf)
   call build_cdf_1d(base_lum_rings, NR, base_ring_cdf)

   p_base = L_base / L_total
   packet_luminosity = L_total / real(N_PACKETS,dp)

   u_thr = 0.0d0
   Fr_thr = 0.0d0
   Fz_thr = 0.0d0
   visit_thr = 0.0d0
   side_thr = 0.0d0
   scatter_thr = 0

   !$omp parallel default(shared) private(k)
   call mc_parallel_loop()
   !$omp end parallel

   u_accum = sum(u_thr, dim=3)
   Fr_accum = sum(Fr_thr, dim=3)
   Fz_accum = sum(Fz_thr, dim=3)

   do i = 1, NR
      do j = 1, NZ
         u_rad(i,j) = u_accum(i,j) / cell_vol(i,j)
         Fr(i,j)    = Fr_accum(i,j) / cell_vol(i,j)
         Fz(i,j)    = Fz_accum(i,j) / cell_vol(i,j)
      end do
   end do
   Fr(1,:) = 0.0d0

   ! ========================================================
   ! Temporal averaging of radiation energy density
   ! This is the field used for pressure and velocity correction.
   ! ========================================================
   if (.not. have_u_rad_avg) then
      u_rad_avg = u_rad
      have_u_rad_avg = .true.
   else
      u_rad_avg = (1.0d0 - TEMP_AVG_ALPHA) * u_rad_avg &
                + TEMP_AVG_ALPHA * u_rad
   end if
   u_rad_used = u_rad_avg

   do i = 1, NR
      do j = 1, NZ
         T_K(i,j)   = max((u_rad_used(i,j)/a_rad)**0.25d0, 1.0d0)
         T_keV(i,j) = k_B * T_K(i,j) / erg_keV
         P_gas(i,j) = 2.0d0 * ne_grid(i,j) * k_B * T_K(i,j)
         P_rad(i,j) = u_rad_used(i,j) / 3.0d0
      end do
   end do

   call gradient_z(P_gas, dPgas_dz, dz)
   call gradient_z(P_rad, dPrad_dz, dz)

   do i = 1, NR
      do j = 1, NZ
         dv_dz_suggested(i,j) = - ( g_grid(i,j) + dPgas_dz(i,j)/rho_grid(i,j) &
                                    + dPrad_dz(i,j)/rho_grid(i,j) ) &
                                / max(v_grid(i,j), 1.0d-30)
      end do
   end do

   ! ========================================================
   ! Iteration residual:
   ! current velocity gradient versus target velocity gradient
   ! ========================================================
   do i = 1, NR
      do j = 1, NZ
         dv_dz_residual(i,j) = &
              (dv_dz_grid(i,j) - dv_dz_suggested(i,j)) &
              / max(abs(dv_dz_suggested(i,j)), RES_EPS)
      end do
   end do
   iter_metric = sqrt(sum(dv_dz_residual**2) / real(NR*NZ, dp))
   write(*,'(A,ES12.4)') "Velocity-gradient residual metric = ", iter_metric

   v_target = 0.0d0
   v_target(:,NZ) = v_ff_grid(:,NZ)

   do j = NZ-1, 1, -1
      do i = 1, NR
         v_target(i,j) = v_target(i,j+1) - dv_dz_suggested(i,j+1) * dz
      end do
   end do

   do i = 1, NR
      do j = 1, NZ
         v_target(i,j) = max(v_target(i,j), V_MIN)
         v_target(i,j) = min(v_target(i,j), v_ff_grid(i,j))

         !========================================================
         ! Under-relaxation with limited velocity correction
         !========================================================
         dv_raw = VELOCITY_RELAX * (v_target(i,j) - v_grid(i,j))

         ! Limit change to 10 percent per iteration.
         dv_lim = 0.10d0 * v_grid(i,j)
         if (dv_raw > dv_lim) then
            dv_raw = dv_lim
         end if
         if (dv_raw < -dv_lim) then
            dv_raw = -dv_lim
         end if
         v_new(i,j) = v_grid(i,j) + dv_raw
      end do
   end do

   call smooth_z(v_new)

   v_grid = v_new
   call write_velocity_profile(iteration)
   scatter_hist = scatter_hist + sum(scatter_thr, dim=2)
end do

call write_outputs()
deallocate(u_thr, Fr_thr, Fz_thr, visit_thr, side_thr, scatter_thr)

contains

  !==================================================================
  !==================================================================
  subroutine write_velocity_profile(iter)
    integer, intent(in) :: iter
    integer :: jj, unit
    character(len=256) :: fname

    write(fname,'("./res/AC_gif/velocity_profile_iter_",I4.4,".dat")') iter
    open(newunit=unit, file=fname, status="replace", action="write")
    write(unit,'(A)') "# z  v_center  v_r_one_third  v_freefall  residual_r_one_third"
    do jj = 1, NZ
       write(unit,'(5ES24.15)') z_cent(jj), &
                              v_grid(1,jj), &
                              v_grid(i_vprof,jj), &
                              v_ff_profile(jj), &
                              dv_dz_residual(i_vprof,jj)
    end do
    close(unit)
  end subroutine write_velocity_profile

  real(dp) function gravity(z)
    real(dp), intent(in) :: z
    gravity = G * M_NS / (R_NS + z)**2
  end function gravity

  real(dp) function v_freefall(z)
    real(dp), intent(in) :: z
    v_freefall = sqrt(2.0d0 * G * M_NS / (R_NS + z))
  end function v_freefall

  real(dp) function rng_uniform(state)
    integer(int64), intent(inout) :: state
    integer(int64) :: x

    x = state
    x = ieor(x, ishft(x, 13))
    x = ieor(x, ishft(x, -7))
    x = ieor(x, ishft(x, 17))
    state = x
    rng_uniform = real(iand(x, huge(1_int64)), dp) / real(huge(1_int64), dp)
    if (rng_uniform <= 0.0d0) rng_uniform = 1.0d-16
    if (rng_uniform >= 1.0d0) rng_uniform = 1.0d0 - 1.0d-16
  end function rng_uniform

  subroutine isotropic_direction_rng(state, dir)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: dir(3)
    real(dp) :: mu, phi, st

    mu = 2.0d0*rng_uniform(state) - 1.0d0
    phi = 2.0d0*pi*rng_uniform(state)
    st = sqrt(max(0.0d0, 1.0d0 - mu*mu))
    dir(1) = st*cos(phi)
    dir(2) = st*sin(phi)
    dir(3) = mu
  end subroutine isotropic_direction_rng

  subroutine upward_direction_rng(state, dir)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: dir(3)
    real(dp) :: mu, phi, st

    mu = rng_uniform(state)
    phi = 2.0d0*pi*rng_uniform(state)
    st = sqrt(max(0.0d0, 1.0d0 - mu*mu))
    dir(1) = st*cos(phi)
    dir(2) = st*sin(phi)
    dir(3) = mu
  end subroutine upward_direction_rng

  subroutine lambert_upward_direction_rng(state, dir)
  integer(int64), intent(inout) :: state
  real(dp), intent(out) :: dir(3)
  real(dp) :: mu, phi, st
    ! Lambert law:
    ! p(mu) = 2 mu, 0 <= mu <= 1
    ! CDF(mu) = mu^2
    ! Therefore mu = sqrt(xi).
    mu = sqrt(rng_uniform(state))
    phi = 2.0d0*pi*rng_uniform(state)
    st = sqrt(max(0.0d0, 1.0d0 - mu*mu))
    dir(1) = st*cos(phi)
    dir(2) = st*sin(phi)
    dir(3) = mu
  end subroutine lambert_upward_direction_rng

  subroutine build_cdf_2d(arr, n1, n2, cdf)
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: arr(n1,n2)
    real(dp), intent(out) :: cdf(n1*n2)
    integer :: ii, jj, kk
    real(dp) :: total

    total = sum(arr)
    kk = 0
    do ii = 1, n1
       do jj = 1, n2
          kk = kk + 1
          if (kk == 1) then
             cdf(kk) = arr(ii,jj) / total
          else
             cdf(kk) = cdf(kk-1) + arr(ii,jj) / total
          end if
       end do
    end do
    cdf(n1*n2) = 1.0d0
  end subroutine build_cdf_2d

  subroutine build_cdf_1d(arr, n, cdf)
    integer, intent(in) :: n
    real(dp), intent(in) :: arr(n)
    real(dp), intent(out) :: cdf(n)
    integer :: ii
    real(dp) :: total

    total = sum(arr)
    do ii = 1, n
       if (ii == 1) then
          cdf(ii) = arr(ii) / total
       else
          cdf(ii) = cdf(ii-1) + arr(ii) / total
       end if
    end do
    cdf(n) = 1.0d0
  end subroutine build_cdf_1d

  integer function sample_cdf(cdf, n, u)
    integer, intent(in) :: n
    real(dp), intent(in) :: cdf(n), u
    integer :: lo, hi, mid

    lo = 1
    hi = n
    do while (lo < hi)
       mid = (lo + hi) / 2
       if (u <= cdf(mid)) then
          hi = mid
       else
          lo = mid + 1
       end if
    end do
    sample_cdf = lo
  end function sample_cdf

  subroutine sample_volume_source(state, pos)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: pos(3)
    integer :: kk, ii, jj
    real(dp) :: r1, r2, rr, phi

    kk = sample_cdf(cell_cdf, NR*NZ, rng_uniform(state))
    ii = (kk-1)/NZ + 1
    jj = mod(kk-1,NZ) + 1

    r1 = r_edges(ii-1)
    r2 = r_edges(ii)

    rr = sqrt(r1*r1 + rng_uniform(state)*(r2*r2 - r1*r1))
    phi = 2.0d0*pi*rng_uniform(state)

    pos(1) = rr*cos(phi)
    pos(2) = rr*sin(phi)
    pos(3) = z_edges(jj-1) + rng_uniform(state)*dz
  end subroutine sample_volume_source

  subroutine sample_base_source(state, pos)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: pos(3)
    integer :: ii
    real(dp) :: r1, r2, rr, phi

    ii = sample_cdf(base_ring_cdf, NR, rng_uniform(state))

    r1 = r_edges(ii-1)
    r2 = r_edges(ii)

    rr = sqrt(r1*r1 + rng_uniform(state)*(r2*r2 - r1*r1))
    phi = 2.0d0*pi*rng_uniform(state)

    pos(1) = rr*cos(phi)
    pos(2) = rr*sin(phi)
    pos(3) = 1.0d-5
  end subroutine sample_base_source

  subroutine get_cell_indices(pos, ii, jj, inside)
    real(dp), intent(in) :: pos(3)
    integer, intent(out) :: ii, jj
    logical, intent(out) :: inside
    real(dp) :: rr

    rr = sqrt(pos(1)**2 + pos(2)**2)

    if (rr < 0.0d0 .or. rr >= R_COL .or. pos(3) < 0.0d0 .or. pos(3) >= H_COL) then
       inside = .false.
       ii = -1
       jj = -1
       return
    end if

    ii = min(int(rr/dr) + 1, NR)
    jj = min(int(pos(3)/dz) + 1, NZ)
    inside = .true.
  end subroutine get_cell_indices

  real(dp) function distance_to_z_boundary(pos, dir)
    real(dp), intent(in) :: pos(3), dir(3)

    if (dir(3) > 0.0d0) then
       distance_to_z_boundary = (H_COL - pos(3)) / dir(3)
    else if (dir(3) < 0.0d0) then
       distance_to_z_boundary = -pos(3) / dir(3)
    else
       distance_to_z_boundary = huge(1.0d0)
    end if
  end function distance_to_z_boundary

  real(dp) function distance_to_outer_cylinder(pos, dir)
    real(dp), intent(in) :: pos(3), dir(3)
    real(dp) :: A, B, Cq, D, sqrtD, s1, s2, smin

    A = dir(1)**2 + dir(2)**2
    if (A < 1.0d-14) then
       distance_to_outer_cylinder = huge(1.0d0)
       return
    end if

    B = 2.0d0 * (pos(1)*dir(1) + pos(2)*dir(2))
    Cq = pos(1)**2 + pos(2)**2 - R_COL**2
    D = B**2 - 4.0d0*A*Cq

    if (D <= 0.0d0) then
       distance_to_outer_cylinder = huge(1.0d0)
       return
    end if

    sqrtD = sqrt(D)
    s1 = (-B + sqrtD) / (2.0d0*A)
    s2 = (-B - sqrtD) / (2.0d0*A)

    smin = huge(1.0d0)
    if (s1 > 1.0d-10) smin = min(smin, s1)
    if (s2 > 1.0d-10) smin = min(smin, s2)

    distance_to_outer_cylinder = smin
  end function distance_to_outer_cylinder

  subroutine gradient_z(arr, grad, dzloc)
    real(dp), intent(in) :: arr(NR,NZ), dzloc
    real(dp), intent(out) :: grad(NR,NZ)
    integer :: ii, jj

    do ii = 1, NR
       grad(ii,1) = (arr(ii,2) - arr(ii,1)) / dzloc
       do jj = 2, NZ-1
          grad(ii,jj) = (arr(ii,jj+1) - arr(ii,jj-1)) / (2.0d0*dzloc)
       end do
       grad(ii,NZ) = (arr(ii,NZ) - arr(ii,NZ-1)) / dzloc
    end do
  end subroutine gradient_z

  subroutine smooth_z(arr)
    real(dp), intent(inout) :: arr(NR,NZ)
    real(dp) :: tmp(NR,NZ)
    integer :: ii, jj

    tmp = arr
    do ii = 1, NR
       do jj = 2, NZ-1
          arr(ii,jj) = 0.25d0*tmp(ii,jj-1) + 0.5d0*tmp(ii,jj) + 0.25d0*tmp(ii,jj+1)
       end do
    end do
  end subroutine smooth_z

  !==================================================================
  ! Relativistic aberration of a photon direction.
  !
  ! n_in    : photon direction in the original frame.
  ! beta_vec: velocity of the new frame relative to the original frame,
  !           in units of c.
  ! n_out   : photon direction in the new frame.
  !
  ! Formula:
  ! n'_parallel = (n_parallel - beta) / (1 - beta dot n)
  ! n'_perp     = n_perp / [gamma * (1 - beta dot n)]
  !
  ! Here beta_vec is the boost from the original frame to the new frame.
  !==================================================================
  subroutine aberrate_direction(n_in, beta_vec, n_out)
    real(dp), intent(in) :: n_in(3)
    real(dp), intent(in) :: beta_vec(3)
    real(dp), intent(out) :: n_out(3)

    real(dp) :: beta2, beta, gamma
    real(dp) :: bdotn, denom
    real(dp) :: npar_scalar
    real(dp) :: npar(3), nperp(3), bhat(3)

    beta2 = dot_product(beta_vec, beta_vec)

    if (beta2 <= 1.0d-30) then
       n_out = n_in
       return
    end if

    beta2 = min(beta2, 0.999999999999d0)
    beta = sqrt(beta2)
    gamma = 1.0d0 / sqrt(1.0d0 - beta2)

    bhat = beta_vec / beta
    bdotn = dot_product(beta_vec, n_in)
    denom = 1.0d0 - bdotn

    if (abs(denom) < 1.0d-30) then
       n_out = n_in
       return
    end if

    npar_scalar = dot_product(n_in, bhat)
    npar = npar_scalar * bhat
    nperp = n_in - npar

    n_out = nperp / (gamma * denom) + &
            ((npar_scalar - beta) / denom) * bhat

    n_out = n_out / sqrt(max(dot_product(n_out, n_out), 1.0d-300))
  end subroutine aberrate_direction

  !==================================================================
  ! Scattering prescription:
  ! 1. Transform incoming lab-frame photon direction to the local
  !    comoving frame of the gas.
  ! 2. Choose the outgoing direction isotropically in the comoving frame.
  ! 3. Transform this outgoing direction back to the lab frame.
  !
  ! The photon energy/packet luminosity is not changed here. This is the
  ! minimal correction needed to make scattering isotropic in the local
  ! zero-velocity frame rather than in the lab frame.
  !==================================================================
  subroutine scatter_isotropic_comoving(state, dir_lab, beta_fluid_lab)
    integer(int64), intent(inout) :: state
    real(dp), intent(inout) :: dir_lab(3)
    real(dp), intent(in) :: beta_fluid_lab(3)

    real(dp) :: dir_com(3), dir_new_com(3), dir_new_lab(3)

    ! Lab frame -> fluid comoving frame.
    call aberrate_direction(dir_lab, beta_fluid_lab, dir_com)

    ! Isotropic Thomson scattering in the comoving frame.
    call isotropic_direction_rng(state, dir_new_com)

    ! Fluid comoving frame -> lab frame.
    call aberrate_direction(dir_new_com, -beta_fluid_lab, dir_new_lab)

    dir_lab = dir_new_lab / sqrt(max(dot_product(dir_new_lab, dir_new_lab), 1.0d-300))
  end subroutine scatter_isotropic_comoving

  subroutine mc_parallel_loop()
    integer :: tid, n, ii, jj, jwall
    integer :: scatter_count, step_count
    integer(int64) :: state
    real(dp) :: pos(3), dir(3), mid(3), er(3)
    real(dp) :: tau_to_scatter, alpha, alpha_lab, s_scatter, s_z, s_side, s_escape, s_next
    real(dp) :: r_mid, r_now, z_now, mu_r, mu_z
    real(dp) :: beta, gamma, beta_vec(3), bdotn
    logical :: inside

    tid = omp_get_thread_num() + 1
    state = seeds(tid)

    !$omp do schedule(dynamic)
    do n = 1, N_PACKETS

       if (rng_uniform(state) < p_base) then
          call sample_base_source(state, pos)
          call upward_direction_rng(state, dir)
       else
          call sample_volume_source(state, pos)
          call isotropic_direction_rng(state, dir)
       end if

       tau_to_scatter = -log(rng_uniform(state))

       scatter_count = 0
       step_count = 0

       do

          step_count = step_count + 1
          if (step_count > MAX_STEPS) exit

          call get_cell_indices(pos, ii, jj, inside)
          if (.not. inside) exit

          visit_thr(ii,jj,tid) = visit_thr(ii,jj,tid) + 1.0d0

          alpha = alpha_grid(ii,jj)

          ! ==========================================================
          ! Relativistic correction to the scattering mean free path.
          !
          ! alpha_grid is interpreted as the Thomson scattering
          ! coefficient in the local comoving frame of the gas.
          !
          ! The gas falls downward. Since +z points upward from the
          ! neutron-star surface, the lab-frame fluid velocity is -z.
          ! v_grid is treated as a positive speed magnitude.
          !
          ! The lab-frame scattering coefficient along photon direction
          ! n is alpha_lab = alpha_com * gamma * (1 - beta dot n).
          ! ==========================================================
          beta = min(v_grid(ii,jj) / c_light, 0.999999d0)
          gamma = 1.0d0 / sqrt(1.0d0 - beta*beta)

          beta_vec(1) = 0.0d0
          beta_vec(2) = 0.0d0
          beta_vec(3) = -beta

          bdotn = dot_product(beta_vec, dir)
          alpha_lab = alpha * gamma * (1.0d0 - bdotn)

          if (alpha_lab > 0.0d0) then
             s_scatter = tau_to_scatter / alpha_lab
          else
             s_scatter = huge(1.0d0)
          end if

          s_z    = distance_to_z_boundary(pos, dir)
          s_side = distance_to_outer_cylinder(pos, dir)

          s_escape = min(s_z, s_side)
          s_next   = min(s_scatter, s_escape)

          if (.not. ieee_is_finite_local(s_next)) exit

          if (s_next < 1.0d-12) then
             pos = pos + 1.0d-5 * dir
             cycle
          end if

          mid = pos + 0.5d0*s_next*dir
          r_mid = sqrt(mid(1)**2 + mid(2)**2)

          if (r_mid > 0.5d0*dr) then
             er(1) = mid(1) / r_mid
             er(2) = mid(2) / r_mid
             er(3) = 0.0d0
             mu_r = dir(1)*er(1) + dir(2)*er(2)
          else
             mu_r = 0.0d0
          end if

          mu_z = dir(3)

          u_thr(ii,jj,tid)  = u_thr(ii,jj,tid)  + packet_luminosity*s_next/c_light
          Fr_thr(ii,jj,tid) = Fr_thr(ii,jj,tid) + packet_luminosity*mu_r*s_next
          Fz_thr(ii,jj,tid) = Fz_thr(ii,jj,tid) + packet_luminosity*mu_z*s_next

          pos = pos + s_next*dir

          if (s_escape <= s_scatter) then

             r_now = sqrt(pos(1)**2 + pos(2)**2)
             z_now = pos(3)

             if (r_now >= R_COL*(1.0d0 - 1.0d-8)) then

                if (z_now >= 0.0d0 .and. z_now < H_COL) then
                   jwall = min(int(z_now/dz) + 1, NZ)
                   side_thr(jwall,tid) = side_thr(jwall,tid) + packet_luminosity
                end if

                exit

             else if (z_now >= H_COL*(1.0d0 - 1.0d-8)) then
                exit
             else if (z_now <= 1.0d-10) then
               call lambert_upward_direction_rng(state, dir)
               pos(3) = 1.0d-5
               cycle
             end if

          end if

          scatter_count = scatter_count + 1

          if (scatter_count > MAX_SCATTERS) exit

          ! ==========================================================
          ! Old version:
          !   call isotropic_direction_rng(state, dir)
          !
          ! New version:
          !   scattering is isotropic only in the local comoving frame.
          !   In the lab frame this produces the correct SR anisotropy.
          ! ==========================================================
          call scatter_isotropic_comoving(state, dir, beta_vec)
          tau_to_scatter = -log(rng_uniform(state))

       end do

       if (scatter_count > MAX_SCATTERS) then
          scatter_thr(MAX_SCATTERS+1,tid) = scatter_thr(MAX_SCATTERS+1,tid) + 1
       else
          scatter_thr(scatter_count,tid) = scatter_thr(scatter_count,tid) + 1
       end if

    end do
    !$omp end do

    seeds(tid) = state

  end subroutine mc_parallel_loop

  logical function ieee_is_finite_local(x)
    real(dp), intent(in) :: x
    ieee_is_finite_local = (abs(x) < huge(1.0d0))
  end function ieee_is_finite_local

  subroutine write_outputs()
    integer :: ii, jj, unit
    real(dp) :: F_side_esc_z(NZ)
    real(dp) :: visit_final(NR,NZ)

    F_side_esc_z = sum(side_thr, dim=2) / side_area_zbin
    visit_final = sum(visit_thr, dim=3)

    open(newunit=unit, file="./res/text_AC", status="replace", action="write")
    do jj = 1, NZ
       write(unit,'(6ES24.15)') z_cent(jj), v_grid(1,jj), v_grid(i_vprof,jj), &
                              dv_dz_residual(i_vprof,jj), &
                              dv_dz_grid(i_vprof,jj), dv_dz_suggested(i_vprof,jj)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC2d", status="replace", action="write")
    do ii = 1, NR
       do jj = 1, NZ
          write(unit,'(3ES24.15)') r_cent(ii), z_cent(jj), v_grid(ii,jj)
       end do
       write(unit,*)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_maps", status="replace", action="write")
    do ii = 1, NR
       do jj = 1, NZ
          write(unit,'(9ES24.15)') r_cent(ii), z_cent(jj), q_grid(ii,jj), &
               u_rad(ii,jj), Fr(ii,jj), Fz(ii,jj), T_keV(ii,jj), &
               P_gas(ii,jj), P_rad(ii,jj)
       end do
       write(unit,*)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_side_flux", status="replace", action="write")
    do jj = 1, NZ
       write(unit,'(2ES24.15)') z_cent(jj), F_side_esc_z(jj)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_scatter_hist", status="replace", action="write")
    do ii = 0, MAX_SCATTERS+1
       write(unit,'(I10,1X,I16)') ii, scatter_hist(ii)
    end do
    close(unit)

    open(newunit=unit, file="./res/test_AC_visit", status="replace", action="write")
    do ii = 1, NR
       do jj = 1, NZ
          write(unit,'(3ES24.15)') r_cent(ii), z_cent(jj), visit_final(ii,jj)
       end do
       write(unit,*)
    end do
    close(unit)

    write(*,*)
    write(*,*) "Output written to:"
    write(*,*) "./res/text_AC"
    write(*,*) "./res/test_AC2d"
    write(*,*) "./res/test_AC_maps"
    write(*,*) "./res/test_AC_side_flux"
    write(*,*) "./res/test_AC_scatter_hist"
    write(*,*) "./res/test_AC_visit"

  end subroutine write_outputs

end subroutine test_AC



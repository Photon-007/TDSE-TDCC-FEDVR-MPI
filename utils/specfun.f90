MODULE specfun
  USE utils
  USE nrtype, only : PI
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dpv = SELECTED_REAL_KIND(12, 60)
CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE laguerre(nd,nu,x,lag)
 IMPLICIT NONE
 INTEGER,INTENT(IN) :: nd, nu
 REAL(dpv), INTENT(OUT) :: lag
 REAL(dpv), INTENT(IN OUT) :: x
 REAL(dpv) :: hg, a, b, posh
 a = -1.d0 * nd
 b = 1.d0 + nu
 lag = 0.0d0
 call chgm(a,b,x,hg)
 call factorial(nd,a)
 call poshsym(B,nd,posh)
 lag = hg * posh / a
END SUBROUTINE laguerre
!------------------------------------------------------------------------------
SUBROUTINE spharm(l,m,theta,phi,cback)
IMPLICIT NONE
 INTEGER, INTENT(IN) :: l, m
 REAL(dpv), INTENT(IN) :: phi, theta
 REAL(dpv) :: x, fact1, fact2
 COMPLEX(dpv), INTENT(OUT) :: cback
 REAL(dpv), DIMENSION(:,:), ALLOCATABLE :: pm,pd
 INTEGER :: mm
 mm = MAX(l,abs(m)) + 1
 x = dcos(theta)
 call factorial(l-abs(m),fact1)
 call factorial(l+abs(m),fact2)
 allocate(pm(0:mm,0:mm),pd(0:mm,0:mm))
 call lpmn(mm, abs(m), l, x, pm, pd)
 cback = cdexp(IU*m*phi)*dsqrt(fact1*(2*l+1)/fact2/4.d0/PI)*pm(abs(m),l)
 if(m .LT. 0)then
    cback = ((-1)**m)*cback
 endif
END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE chgm(a, b, x, hg)

!       ===================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Routine called: GAMMA for computing â(x)
!       ===================================================

REAL (dpv), INTENT(IN OUT)  :: a
REAL (dpv), INTENT(IN)      :: b
REAL (dpv), INTENT(IN OUT)  :: x
REAL (dpv), INTENT(OUT)     :: hg

REAL (dpv), PARAMETER  :: pi = 3.141592653589793_dpv
REAL (dpv)  :: a0, a1, hg1, hg2, r, r1, r2, rg, sum1, sum2, ta, tb, tba,  &
              x0, xg, y0, y1
INTEGER    :: i, j, k, la, m, n, nl

a0 = a
a1 = a
x0 = x
hg = 0.0_dpv
IF (b == 0.0_dpv .OR. b == -ABS(INT(b))) THEN
  hg = 1.0D+300
ELSE IF (a == 0.0_dpv .OR. x == 0.0_dpv) THEN
  hg = 1.0_dpv
ELSE IF (a == -1.0_dpv) THEN
  hg = 1.0_dpv - x / b
ELSE IF (a == b) THEN
  hg = EXP(x)
ELSE IF (a-b == 1.0_dpv) THEN
  hg = (1.0_dpv + x/b) * EXP(x)
ELSE IF (a == 1.0_dpv .AND. b == 2.0_dpv) THEN
  hg = (EXP(x)-1.0_dpv) / x
ELSE IF (a == INT(a) .AND. a < 0.0_dpv) THEN
  m = INT(-a)
  r = 1.0_dpv
  hg = 1.0_dpv
  DO  k = 1, m
    r = r * (a+k-1) / k / (b+k-1) * x
    hg = hg + r
  END DO
END IF
IF (hg /= 0.0_dpv) RETURN
IF (x < 0.0_dpv) THEN
  a = b - a
  a0 = a
  x = ABS(x)
END IF
IF (a < 2.0_dpv) nl = 0
IF (a >= 2.0_dpv) THEN
  nl = 1
  la = INT(a)
  a = a - la - 1
END IF
DO  n = 0, nl
  IF (a0 >= 2.0_dpv) a = a + 1.0_dpv
  IF (x <= 30.0_dpv + ABS(b) .OR. a < 0.0_dpv) THEN
    hg = 1.0_dpv
    rg = 1.0_dpv
    DO  j = 1, 500
      rg = rg * (a+j-1) / (j*(b+j-1)) * x
      hg = hg + rg
      IF (ABS(rg/hg) < 1.0D-15) GO TO 40
    END DO
  ELSE
    CALL gamma(a, ta)
    CALL gamma(b, tb)
    xg = b - a
    CALL gamma(xg, tba)
    sum1 = 1.0_dpv
    sum2 = 1.0_dpv
    r1 = 1.0_dpv
    r2 = 1.0_dpv
    DO  i = 1, 8
      r1 = -r1 * (a+i-1) * (a-b+i) / (x*i)
      r2 = -r2 * (b-a+i-1) * (a-i) / (x*i)
      sum1 = sum1 + r1
      sum2 = sum2 + r2
    END DO
    hg1 = tb / tba * x ** (-a) * COS(pi*a) * sum1
    hg2 = tb / ta * EXP(x) * x ** (a-b) * sum2
    hg = hg1 + hg2
  END IF
  40 IF (n == 0) y0 = hg
  IF (n == 1) y1 = hg
END DO
IF (a0 >= 2.0_dpv) THEN
  DO  i = 1, la - 1
    hg = ((2.0_dpv*a-b+x)*y1 + (b-a)*y0) / a
    y0 = y1
    y1 = hg
    a = a + 1.0_dpv
  END DO
END IF
IF (x0 < 0.0_dpv) hg = hg * EXP(x0)
a = a1
x = x0
RETURN
END SUBROUTINE chgm

SUBROUTINE besselj(l,x,bjval,byval)
  IMPLICIT NONE
  REAL(dpv), INTENT(IN) :: x
  INTEGER, INTENT(IN) :: l
  REAL(dpv), INTENT(OUT) :: bjval, byval
  REAL(dpv), DIMENSION(:), ALLOCATABLE :: bj, dj, by, dy
  REAL(DPV) :: vm
  REAL(dpv) :: v
  v = 0.5d0 + l

  if(allocated(bj))deallocate(bj)
  if(allocated(dj))deallocate(dj)
  if(allocated(by))deallocate(by)
  if(allocated(dy))deallocate(dy)

  allocate(bj(0:l),dj(0:l),by(0:l),dy(0:l))

  ! write(6,*)"CP2"
  call jyv(v, x, vm, bj, dj, by, dy)
  ! write(6,*)x,bj(l)
  bjval = bj(l)*dsqrt(PI/2.d0/x)
  byval = by(l)*dsqrt(PI/2.d0/x)
  deallocate(bj,dj,by,dy)
END SUBROUTINE besselj

SUBROUTINE jyv(v, x, vm, bj, dj, by, dy)

!    =======================================================
!    Purpose: Compute Bessel functions Jv(x) and Yv(x)
!             and their derivatives
!    Input :  x --- Argument of Jv(x) and Yv(x)
!             v --- Order of Jv(x) and Yv(x)
!                   ( v = n+v0, 0 ó v0 < 1, n = 0,1,2,... )
!    Output:  BJ(n) --- Jn+v0(x)
!             DJ(n) --- Jn+v0'(x)
!             BY(n) --- Yn+v0(x)
!             DY(n) --- Yn+v0'(x)
!             VM --- Highest order computed
!    Routines called:
!         (1) GAMMA for computing gamma function
!         (2) MSTA1 and MSTA2 for computing the starting
!             point for backward recurrence
!    =======================================================

REAL (dpv), INTENT(IN)   :: v
REAL (dpv), INTENT(IN)   :: x
REAL (dpv), INTENT(OUT)  :: vm
REAL (dpv), INTENT(OUT)  :: bj(0:)
REAL (dpv), INTENT(OUT)  :: dj(0:)
REAL (dpv), INTENT(OUT)  :: by(0:)
REAL (dpv), INTENT(OUT)  :: dy(0:)

REAL (dpv), PARAMETER  :: el = .5772156649015329_dpv, pi = 3.141592653589793_dpv, &
                         rp2 = .63661977236758_dpv
REAL (dpv)  :: a, a0, b, bju0, bju1, bjv0, bjv1, bjvl, byv0, byv1, byvk,  &
              ck, cs, cs0, cs1, ec, f, f0, f1, f2, ga, gb, pv0, pv1, px,  &
              qx, r, r0, r1, rp, rq, sk, v0, vg, vl, vv, w0, w1, x2, xk
INTEGER    :: j, k, k0, l, m, n


x2 = x * x
n = v
v0 = v - n
IF (x < 1.0D-100) THEN
  DO  k = 0, n
    bj(k) = 0.0D0
    dj(k) = 0.0D0
    by(k) = -1.0D+300
    dy(k) = 1.0D+300
  END DO
  IF (v0 == 0.0) THEN
    bj(0) = 1.0D0
    dj(1) = 0.5D0
  ELSE
    dj(0) = 1.0D+300
  END IF
  vm = v
  RETURN
END IF

IF (x <= 12.0) THEN
  DO  l = 0, 1
    vl = v0 + l
    bjvl = 1.0D0
    r = 1.0D0
    DO  k = 1, 40
      r = -0.25D0 * r * x2 / (k*(k+vl))
      bjvl = bjvl + r
      IF (ABS(r) < ABS(bjvl)*1.0D-15) EXIT
    END DO
    vg = 1.0D0 + vl
    CALL gamma(vg, ga)
    a = (0.5D0*x) ** vl / ga
    IF (l == 0) bjv0 = bjvl * a
    IF (l == 1) bjv1 = bjvl * a
  END DO
ELSE
  k0 = 11
  IF (x >= 35.0) k0 = 10
  IF (x >= 50.0) k0 = 8
  DO  j = 0, 1
    vv = 4.0D0 * (j+v0) * (j+v0)
    px = 1.0D0
    rp = 1.0D0
    DO  k = 1, k0
      rp = -0.78125D-2 * rp * (vv - (4*k-3)**2) * (vv - (4*k-1)**2) /  &
           (k*(2*k-1)*x2)
      px = px + rp
    END DO
    qx = 1.0D0
    rq = 1.0D0
    DO  k = 1, k0
      rq = -0.78125D-2 * rq * (vv - (4*k-1)**2) * (vv - (4*k+1)**2) /   &
           (k*(2*k+1)*x2)
      qx = qx + rq
    END DO
    qx = 0.125D0 * (vv-1.0) * qx / x
    xk = x - (0.5D0*(j+v0) + 0.25D0) * pi
    a0 = SQRT(rp2/x)
    ck = COS(xk)
    sk = SIN(xk)
    IF (j == 0) THEN
      bjv0 = a0 * (px*ck - qx*sk)
      byv0 = a0 * (px*sk + qx*ck)
    ELSE IF (j == 1) THEN
      bjv1 = a0 * (px*ck - qx*sk)
      byv1 = a0 * (px*sk + qx*ck)
    END IF
  END DO
END IF


bj(0) = bjv0
bj(1) = bjv1
dj(0) = v0 / x * bj(0) - bj(1)
dj(1) = -(1.0D0+v0) / x * bj(1) + bj(0)
IF (n >= 2 .AND. n <= INT(0.9*x)) THEN
 ! write(6,*)x,'cp0a'
  f0 = bjv0
  f1 = bjv1
  DO  k = 2, n
    f = 2.0D0 * (k+v0-1.0D0) / x * f1 - f0
    bj(k) = f
    f0 = f1
    f1 = f
  END DO
ELSE IF (n >= 2) THEN
!  write(6,*)x,'cp0b0'
  m = msta1(x, 200)
  IF (m < n) THEN
    n = m
  ELSE
    if(x .GT. 0.1d0 .and. x .LT. 1.d0)then
       m = msta2(x, n, 7) ! modified fron 15 to 7----means loss of accuracy find the bug latter m -s a large negative number if x\in(0.2,1) !!!!!!!!!!!!
    else
       m = msta2(x, n, 15)
    endif
  END IF

  f2 = 0.0D0
  f1 = 1.0D-100

!  write(6,*)x,'cp0b0',m

  DO  k = m, 0, -1
    f = 2.0D0 * (v0+k+1) / x * f1 - f2
    IF (k <= n) bj(k) = f
    f2 = f1
    f1 = f

  END DO

!  write(6,*)x,'cp0b1',f
  IF (ABS(bjv0) > ABS(bjv1)) THEN
    cs = bjv0 / f
  ELSE
    cs = bjv1 / f2
  END IF
  bj(0:n) = cs * bj(0:n)
END IF
!write(6,*)x,'cp1'

DO  k = 2, n
  dj(k) = -(k+v0) / x * bj(k) + bj(k-1)
END DO
IF (x <= 12.0D0) THEN
  IF (v0 /= 0.0) THEN
    DO  l = 0, 1
      vl = v0 + l
      bjvl = 1.0D0
      r = 1.0D0
      DO  k = 1, 40
        r = -0.25D0 * r * x2 / (k*(k-vl))
        bjvl = bjvl + r
        IF (ABS(r) < ABS(bjvl)*1.0D-15) EXIT
      END DO
      vg = 1.0D0 - vl
      CALL gamma(vg,gb)
      b = (2.0D0/x) ** vl / gb
      IF (l == 0) bju0 = bjvl * b
      IF (l == 1) bju1 = bjvl * b
    END DO
    pv0 = pi * v0
    pv1 = pi * (1.0D0+v0)
    byv0 = (bjv0*COS(pv0) - bju0) / SIN(pv0)
    byv1 = (bjv1*COS(pv1) - bju1) / SIN(pv1)
  ELSE
    ec = LOG(x/2.0D0) + el
    cs0 = 0.0D0
    w0 = 0.0D0
    r0 = 1.0D0
    DO  k = 1, 30
      w0 = w0 + 1.0D0 / k
      r0 = -0.25D0 * r0 / (k*k) * x2
      cs0 = cs0 + r0 * w0
    END DO
    byv0 = rp2 * (ec*bjv0 - cs0)
    cs1 = 1.0D0
    w1 = 0.0D0
    r1 = 1.0D0
    DO  k = 1, 30
      w1 = w1 + 1.0D0 / k
      r1 = -0.25D0 * r1 / (k*(k+1)) * x2
      cs1 = cs1 + r1 * (2.0D0*w1+1.0D0/(k+1.0D0))
    END DO
    byv1 = rp2 * (ec*bjv1 - 1.0D0/x - 0.25D0*x*cs1)
  END IF
END IF
by(0) = byv0
by(1) = byv1
DO  k = 2, n
  byvk = 2.0D0 * (v0+k-1) / x * byv1 - byv0
  by(k) = byvk
  byv0 = byv1
  byv1 = byvk
END DO
dy(0) = v0 / x * by(0) - by(1)
DO  k = 1, n
  dy(k) = -(k+v0) / x * by(k) + by(k-1)
END DO
vm = n + v0
RETURN
END SUBROUTINE jyv


SUBROUTINE gamma(x,ga)

!       ==================================================
!       Purpose: Compute gamma function â(x)
!       Input :  x  --- Argument of â(x)
!                       ( x is not equal to 0,-1,-2,úúú)
!       Output:  GA --- â(x)
!       ==================================================


REAL (dpv), INTENT(IN)      :: x
REAL (dpv), INTENT(OUT)     :: ga

REAL (dpv), PARAMETER  :: g(26) = (/1.0_dpv, 0.5772156649015329_dpv,  &
      -0.6558780715202538_dpv, -0.420026350340952D-1, 0.1665386113822915_dpv,  &
      -0.421977345555443D-1, -.96219715278770D-2, .72189432466630D-2,  &
      -0.11651675918591D-2, -.2152416741149D-3, .1280502823882D-3,  &
      -0.201348547807D-4, -.12504934821D-5, 0.11330272320D-5,  &
      -0.2056338417D-6, .61160950D-8, .50020075D-8, -.11812746D-8,  &
       0.1043427D-9, .77823D-11, -.36968D-11, .51D-12, -.206D-13,   &
       -.54D-14, .14D-14, .1D-15 /)
REAL (dpv), PARAMETER  :: pi = 3.141592653589793_dpv
REAL (dpv)  :: gr, r, z
INTEGER    :: k, m, m1

IF (x == INT(x)) THEN
  IF (x > 0.0_dpv) THEN
    ga = 1.0_dpv
    m1 = x - 1
    DO  k = 2, m1
      ga = ga * k
    END DO
  ELSE
    ga = 1.0D+300
  END IF
ELSE
  IF (ABS(x) > 1.0_dpv) THEN
    z = ABS(x)
    m = INT(z)
    r = 1.0_dpv
    DO  k = 1, m
      r = r * (z-k)
    END DO
    z = z - m
  ELSE
    z = x
  END IF

  gr = g(26)
  DO  k = 25, 1, -1
    gr = gr * z + g(k)
  END DO
  ga = 1.0_dpv / (gr*z)
  IF (ABS(x) > 1.0_dpv) THEN
    ga = ga * r
    IF (x < 0.0_dpv) ga = -pi / (x*ga*SIN(pi*x))
  END IF
END IF
RETURN
END SUBROUTINE gamma

SUBROUTINE lpn(n, x, pn, pd)

!    ===============================================
!    Purpose: Compute Legendre polynomials Pn(x)
!             and their derivatives Pn'(x)
!    Input :  x --- Argument of Pn(x)
!             n --- Degree of Pn(x) ( n = 0,1,...)
!    Output:  PN(n) --- Pn(x)
!             PD(n) --- Pn'(x)
!    ===============================================

INTEGER, INTENT(IN)     :: n
REAL (dpv), INTENT(IN)   :: x
REAL (dpv), INTENT(OUT)  :: pn(0:n)
REAL (dpv), INTENT(OUT)  :: pd(0:n)

REAL (dpv)  :: p0, p1, pf
INTEGER    :: k

pn(0) = 1.0_dpv
pn(1) = x
pd(0) = 0.0_dpv
pd(1) = 1.0_dpv
p0 = 1.0_dpv
p1 = x
DO  k = 2, n
  pf = (2.0_dpv*k-1.0_dpv) / k * x * p1 - (k-1.0_dpv) / k * p0
  pn(k) = pf
  IF (ABS(x) == 1.0_dpv) THEN
    pd(k) = 0.5_dpv * x ** (k+1) * k * (k+1)
  ELSE
    pd(k) = k * (p1 - x*pf) / (1.0_dpv - x*x)
  END IF
  p0 = p1
  p1 = pf
END DO
RETURN
END SUBROUTINE lpn

SUBROUTINE cgama(x, y, kf, gr, gi)

!    ======================================================
!    Purpose: Compute the gamma function â(z) or ln[â(z)]
!             for a complex argument
!    Input :  x  --- Real part of z
!             y  --- Imaginary part of z
!             KF --- Function code
!                    KF=0 for ln[â(z)]
!                    KF=1 for â(z)
!    Output:  GR --- Real part of ln[â(z)] or â(z)
!             GI --- Imaginary part of ln[â(z)] or â(z)
!    ======================================================

REAL (dpv), INTENT(IN OUT)  :: x
REAL (dpv), INTENT(IN OUT)  :: y
INTEGER, INTENT(IN)        :: kf
REAL (dpv), INTENT(OUT)     :: gr
REAL (dpv), INTENT(OUT)     :: gi

REAL (dpv), PARAMETER  :: a(10) = (/ 8.333333333333333D-02,  &
   -2.777777777777778D-03,  7.936507936507937D-04, -5.952380952380952D-04,  &
    8.417508417508418D-04, -1.917526917526918D-03,  6.410256410256410D-03,  &
   -2.955065359477124D-02,  1.796443723688307D-01, -1.39243221690590_dpv /)
REAL (dpv), PARAMETER  :: pi = 3.141592653589793_dpv
REAL (dpv)  :: g0, gi1, gr1, sr, si, t, th, th1, th2, x0, x1, y1, z1, z2
INTEGER    :: j, k, na

x1 = x
y1 = y
IF (y == 0.0_dpv .AND. x == INT(x) .AND. x <= 0.0_dpv) THEN
  gr = 1.0D+300
  gi = 0.0_dpv
  RETURN
ELSE IF (x < 0.0_dpv) THEN
  x1 = x
  y1 = y
  x = -x
  y = -y
END IF
x0 = x
IF (x <= 7.0) THEN
  na = INT(7-x)
  x0 = x + na
END IF
z1 = SQRT(x0*x0 + y*y)
th = ATAN(y/x0)
gr = (x0-.5_dpv) * LOG(z1) - th * y - x0 + 0.5_dpv * LOG(2.0_dpv*pi)
gi = th * (x0 - 0.5_dpv) + y * LOG(z1) - y
DO  k = 1, 10
  t = z1 ** (1-2*k)
  gr = gr + a(k) * t * COS((2*k-1)*th)
  gi = gi - a(k) * t * SIN((2*k-1)*th)
END DO
IF (x <= 7.0) THEN
  gr1 = 0.0_dpv
  gi1 = 0.0_dpv
  DO  j = 0, na - 1
    gr1 = gr1 + .5_dpv * LOG((x+j)**2 + y*y)
    gi1 = gi1 + ATAN(y/(x+j))
  END DO
  gr = gr - gr1
  gi = gi - gi1
END IF
IF (x1 < 0.0_dpv) THEN
  z1 = SQRT(x*x + y*y)
  th1 = ATAN(y/x)
  sr = -SIN(pi*x) * COSH(pi*y)
  si = -COS(pi*x) * SINH(pi*y)
  z2 = SQRT(sr*sr + si*si)
  th2 = ATAN(si/sr)
  IF (sr < 0.0_dpv) th2 = pi + th2
  gr = LOG(pi/(z1*z2)) - gr
  gi = -th1 - th2 - gi
  x = x1
  y = y1
END IF
IF (kf == 1) THEN
  g0 = EXP(gr)
  gr = g0 * COS(gi)
  gi = g0 * SIN(gi)
END IF
RETURN
END SUBROUTINE cgama
!------------------------------------------------------------------------------
 SUBROUTINE cjyna(n, z, nm, cbj, cdj, cby, cdy)

!       =======================================================
!       Purpose: Compute Bessel functions Jn(z), Yn(z) and
!                their derivatives for a complex argument
!       Input :  z --- Complex argument of Jn(z) and Yn(z)
!                n --- Order of Jn(z) and Yn(z)
!       Output:  CBJ(n) --- Jn(z)
!                CDJ(n) --- Jn'(z)
!                CBY(n) --- Yn(z)
!                CDY(n) --- Yn'(z)
!                NM --- Highest order computed
!       Routines called:
!            (1) CJY01 to calculate J0(z), J1(z), Y0(z), Y1(z)
!            (2) MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =======================================================

 INTEGER, INTENT(IN)        :: n
 COMPLEX (dpv), INTENT(IN)   :: z
 INTEGER, INTENT(OUT)       :: nm
 COMPLEX (dpv), INTENT(OUT)  :: cbj(0:)
 COMPLEX (dpv), INTENT(OUT)  :: cdj(0:)
 COMPLEX (dpv), INTENT(OUT)  :: cby(0:)
 COMPLEX (dpv), INTENT(OUT)  :: cdy(0:)

 COMPLEX (dpv)  :: cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1,  &
                 cf, cf1, cf2, cg0, cg1, ch0, ch1, ch2, cj0, cj1, cjk,  &
                 cp11, cp12, cp21, cp22, cs, cyk, cyl1, cyl2, cylk
 REAL (dpv)     :: a0, wa, ya0, ya1, yak
 INTEGER       :: k, lb, lb0, m
 REAL (dpv), PARAMETER  :: pi = 3.141592653589793_dpv

 a0 = ABS(z)
 nm = n
 IF (a0 < 1.0D-100) THEN
   DO  k = 0, n
     cbj(k) = (0.0_dpv,0.0_dpv)
     cdj(k) = (0.0_dpv,0.0_dpv)
     cby(k) = -(1.0D+300,0.0_dpv)
     cdy(k) = (1.0D+300,0.0_dpv)
   END DO
   cbj(0) = (1.0_dpv,0.0_dpv)
   cdj(1) = (0.5_dpv,0.0_dpv)
   RETURN
 END IF
 CALL cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)
 cbj(0) = cbj0
 cbj(1) = cbj1
 cby(0) = cby0
 cby(1) = cby1
 cdj(0) = cdj0
 cdj(1) = cdj1
 cdy(0) = cdy0
 cdy(1) = cdy1
 IF (n <= 1) RETURN
 IF (n < INT(0.25*a0)) THEN
   cj0 = cbj0
   cj1 = cbj1
   DO  k = 2, n
     cjk = 2 * (k-1) / z * cj1 - cj0
     cbj(k) = cjk
     cj0 = cj1
     cj1 = cjk
   END DO
 ELSE
   m = msta1(a0, 200)
   IF (m < n) THEN
     nm = m
   ELSE
     m = msta2(a0, n, 15)
   END IF
   cf2 = (0.0_dp,0.0_dpv)
   cf1 = (1.0D-100,0.0_dpv)
   DO  k = m, 0, -1
     cf = 2 * (k+1) / z * cf1 - cf2
     IF (k <= nm) cbj(k) = cf
     cf2 = cf1
     cf1 = cf
   END DO
   IF (ABS(cbj0) > ABS(cbj1)) THEN
     cs = cbj0 / cf
   ELSE
     cs = cbj1 / cf2
   END IF
   DO  k = 0, nm
     cbj(k) = cs * cbj(k)
   END DO
 END IF
 DO  k = 2, nm
   cdj(k) = cbj(k-1) - k / z * cbj(k)
 END DO
 ya0 = ABS(cby0)
 lb = 0
 cg0 = cby0
 cg1 = cby1
 DO  k = 2, nm
   cyk = 2 * (k-1) / z * cg1 - cg0
   IF (ABS(cyk) <= 1.0D+290) THEN
     yak = ABS(cyk)
     ya1 = ABS(cg0)
     IF (yak < ya0 .AND. yak < ya1) lb = k
     cby(k) = cyk
     cg0 = cg1
     cg1 = cyk
   END IF
 END DO
 lb0 = 0
 IF (lb > 4 .AND. AIMAG(z) /= 0.0_dpv) THEN
   70 IF (lb /= lb0) THEN
     ch2 = (1.0_dpv,0.0_dpv)
     ch1 = (0.0_dpv,0.0_dpv)
     lb0 = lb
     DO  k = lb, 1, -1
       ch0 = 2.0_dpv * k / z * ch1 - ch2
       ch2 = ch1
       ch1 = ch0
     END DO
     cp12 = ch0
     cp22 = ch2
     ch2 = (0.0_dpv,0.0_dpv)
     ch1 = (1.0_dpv,0.0_dpv)
     DO  k = lb, 1, -1
       ch0 = 2 * k / z * ch1 - ch2
       ch2 = ch1
       ch1 = ch0
     END DO
     cp11 = ch0
     cp21 = ch2
     IF (lb == nm) cbj(lb+1) = 2.0_dpv * lb / z * cbj(lb) - cbj(lb-1)
     IF (ABS(cbj(0)) > ABS(cbj(1))) THEN
       cby(lb+1) = (cbj(lb+1)*cby0 - 2.0_dpv*cp11/(pi*z)) / cbj(0)
       cby(lb) = (cbj(lb)*cby0 + 2.0_dpv*cp12/(pi*z)) / cbj(0)
     ELSE
       cby(lb+1) = (cbj(lb+1)*cby1 - 2.0_dpv*cp21/(pi*z)) / cbj(1)
       cby(lb) = (cbj(lb)*cby1 + 2.0_dp*cp22/(pi*z)) / cbj(1)
     END IF
     cyl2 = cby(lb+1)
     cyl1 = cby(lb)
     DO  k = lb - 1, 0, -1
       cylk = 2 * (k+1) / z * cyl1 - cyl2
       cby(k) = cylk
       cyl2 = cyl1
       cyl1 = cylk
     END DO
     cyl1 = cby(lb)
     cyl2 = cby(lb+1)
     DO  k = lb + 1, nm - 1
       cylk = 2 * k / z * cyl2 - cyl1
       cby(k+1) = cylk
       cyl1 = cyl2
       cyl2 = cylk
     END DO
     DO  k = 2, nm
       wa = ABS(cby(k))
       IF (wa < ABS(cby(k-1))) lb = k
     END DO
     GO TO 70
   END IF
 END IF
 DO  k = 2, nm
   cdy(k) = cby(k-1) - k / z * cby(k)
 END DO
 RETURN
 END SUBROUTINE cjyna
!------------------------------------------------------------------------------
SUBROUTINE cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)

!       =======================================================
!       Purpose: Compute Bessel functions J0(z), J1(z), Y0(z),
!                Y1(z), and their derivatives for a complex
!                argument
!       Input :  z --- Complex argument
!       Output:  CBJ0 --- J0(z)
!                CDJ0 --- J0'(z)
!                CBJ1 --- J1(z)
!                CDJ1 --- J1'(z)
!                CBY0 --- Y0(z)
!                CDY0 --- Y0'(z)
!                CBY1 --- Y1(z)
!                CDY1 --- Y1'(z)
!       =======================================================

COMPLEX (dpv), INTENT(IN)   :: z
COMPLEX (dpv), INTENT(OUT)  :: cbj0
COMPLEX (dpv), INTENT(OUT)  :: cdj0
COMPLEX (dpv), INTENT(OUT)  :: cbj1
COMPLEX (dpv), INTENT(OUT)  :: cdj1
COMPLEX (dpv), INTENT(OUT)  :: cby0
COMPLEX (dpv), INTENT(OUT)  :: cdy0
COMPLEX (dpv), INTENT(OUT)  :: cby1
COMPLEX (dpv), INTENT(OUT)  :: cdy1

REAL (dpv), PARAMETER  :: a(12) = (/ -0.703125D-01, .112152099609375_dpv,   &
     -.5725014209747314_dpv, .6074042001273483D+01, -  &
      .1100171402692467D+03, .3038090510922384D+04, -  &
      .1188384262567832D+06, .6252951493434797D+07, -  &
      .4259392165047669D+09, .3646840080706556D+11, -  &
      .3833534661393944D+13, .4854014686852901D+15 /)
REAL (dpv), PARAMETER  :: b(12) = (/.732421875D-01, -.2271080017089844_dpv,  &
      .1727727502584457D+01, -.2438052969955606D+02,  &
      .5513358961220206D+03, -.1825775547429318D+05,  &
      .8328593040162893D+06, -.5006958953198893D+08,  &
      .3836255180230433D+10, -.3649010818849833D+12,  &
      .4218971570284096D+14, -.5827244631566907D+16 /)
REAL (dpv), PARAMETER  :: a1(12) = (/ 0.1171875_dpv, -.144195556640625_dpv,  &
      .6765925884246826_dpv, -.6883914268109947D+01,  &
      .1215978918765359D+03, -.3302272294480852D+04,  &
      .1276412726461746D+06, -.6656367718817688D+07,  &
      .4502786003050393D+09, -.3833857520742790D+11,  &
      .4011838599133198D+13, -.5060568503314727D+15 /)
REAL (dpv), PARAMETER  :: b1(12) = (/ -0.1025390625_dpv, .2775764465332031_dpv,  &
     -.1993531733751297D+01, .2724882731126854D+02, -  &
      .6038440767050702D+03, .1971837591223663D+05, -  &
      .8902978767070678D+06, .5310411010968522D+08, -  &
      .4043620325107754D+10, .3827011346598605D+12, -  &
      .4406481417852278D+14, .6065091351222699D+16 /)

REAL (dpv)     :: a0, rp2, w0, w1
INTEGER       :: k, k0
COMPLEX (dpv)  :: ci, cp, cp0, cp1, cq0, cq1, cr, cs, ct1, ct2, cu, z1, z2
REAL (dpv), PARAMETER  :: pi = 3.141592653589793_dpv, el = 0.5772156649015329_dpv

rp2 = 2.0_dpv / pi
ci = (0.0_dpv,1.0_dpv)
a0 = ABS(z)
z2 = z * z
z1 = z
IF (a0 == 0.0_dpv) THEN
  cbj0 = (1.0_dpv,0.0_dpv)
  cbj1 = (0.0_dpv,0.0_dpv)
  cdj0 = (0.0_dpv,0.0_dpv)
  cdj1 = (0.5_dpv,0.0_dpv)
  cby0 = -(1.0D300,0.0_dpv)
  cby1 = -(1.0D300,0.0_dpv)
  cdy0 = (1.0D300,0.0_dpv)
  cdy1 = (1.0D300,0.0_dpv)
  RETURN
END IF
IF (REAL(z) < 0.0) z1 = -z
IF (a0 <= 12.0) THEN
  cbj0 = (1.0_dpv,0.0_dpv)
  cr = (1.0_dpv,0.0_dpv)
  DO  k = 1, 40
    cr = -0.25_dpv * cr * z2 / (k*k)
    cbj0 = cbj0 + cr
    IF (ABS(cr) < ABS(cbj0)*1.0D-15) EXIT
  END DO
  cbj1 = (1.0_dpv,0.0_dpv)
  cr = (1.0_dpv,0.0_dpv)
  DO  k = 1, 40
    cr = -0.25_dpv * cr * z2 / (k*(k+1))
    cbj1 = cbj1 + cr
    IF (ABS(cr) < ABS(cbj1)*1.0D-15) EXIT
  END DO
  cbj1 = 0.5_dpv * z1 * cbj1
  w0 = 0.0_dpv
  cr = (1.0_dpv,0.0_dpv)
  cs = (0.0_dpv,0.0_dpv)
  DO  k = 1, 40
    w0 = w0 + 1.0_dpv / k
    cr = -0.25_dpv * cr / (k*k) * z2
    cp = cr * w0
    cs = cs + cp
    IF (ABS(cp) < ABS(cs)*1.0D-15) EXIT
  END DO
  cby0 = rp2 * (LOG(z1/2.0_dpv)+el) * cbj0 - rp2 * cs
  w1 = 0.0_dpv
  cr = (1.0_dpv,0.0_dpv)
  cs = (1.0_dpv,0.0_dpv)
  DO  k = 1, 40
    w1 = w1 + 1.0_dpv / k
    cr = -0.25_dpv * cr / (k*(k+1)) * z2
    cp = cr * (2.0_dpv*w1+1.0_dpv/(k+1.0_dpv))
    cs = cs + cp
    IF (ABS(cp) < ABS(cs)*1.0D-15) EXIT
  END DO
  cby1 = rp2 * ((LOG(z1/2.0_dpv) + el)*cbj1 - 1.0_dpv/z1 - .25_dpv*z1*cs)
ELSE
  k0 = 12
  IF (a0 >= 35.0) k0 = 10
  IF (a0 >= 50.0) k0 = 8
  ct1 = z1 - .25_dpv * pi
  cp0 = (1.0_dpv,0.0_dpv)
  DO  k = 1, k0
    cp0 = cp0 + a(k) * z1 ** (-2*k)
  END DO
  cq0 = -0.125_dpv / z1
  DO  k = 1, k0
    cq0 = cq0 + b(k) * z1 ** (-2*k-1)
  END DO
  cu = SQRT(rp2/z1)
  cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
  cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
  ct2 = z1 - .75_dpv * pi
  cp1 = (1.0_dpv,0.0_dpv)
  DO  k = 1, k0
    cp1 = cp1 + a1(k) * z1 ** (-2*k)
  END DO
  cq1 = 0.375_dpv / z1
  DO  k = 1, k0
    cq1 = cq1 + b1(k) * z1 ** (-2*k-1)
  END DO
  cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
  cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
END IF
IF (REAL(z) < 0.0) THEN
  IF (AIMAG(z) < 0.0) cby0 = cby0 - 2.0_dpv * ci * cbj0
  IF (AIMAG(z) > 0.0) cby0 = cby0 + 2.0_dpv * ci * cbj0
  IF (AIMAG(z) < 0.0) cby1 = -(cby1 - 2.0_dpv*ci*cbj1)
  IF (AIMAG(z) > 0.0) cby1 = -(cby1 + 2.0_dpv*ci*cbj1)
  cbj1 = -cbj1
END IF
cdj0 = -cbj1
cdj1 = cbj0 - 1.0_dpv / z * cbj1
cdy0 = -cby1
cdy1 = cby0 - 1.0_dpv / z * cby1
RETURN
END SUBROUTINE cjy01
!------------------------------------------------------------------------------
FUNCTION msta1(x, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================

REAL (dpv), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

REAL (dp)  :: a0, f, f0, f1
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
n0 = INT(1.1*a0) + 1
f0 = envj(n0,a0) - mp
n1 = n0 + 5
f1 = envj(n1,a0) - mp
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dpv - f0/f1)
  f = envj(nn,a0) - mp
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn
RETURN
END FUNCTION msta1



FUNCTION msta2(x, n, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================

REAL (dpv), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

REAL (dpv)  :: a0, ejn, f, f0, f1, hmp, obj
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
hmp = 0.5_dpv * mp
ejn = envj(n, a0)
IF (ejn <= hmp) THEN
  obj = mp
  n0 = INT(1.1*a0)
ELSE
  obj = hmp + ejn
  n0 = n
END IF
f0 = envj(n0,a0) - obj
n1 = n0 + 5
f1 = envj(n1,a0) - obj
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dpv - f0/f1)
  f = envj(nn, a0) - obj
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn + 10
RETURN
END FUNCTION msta2



FUNCTION envj(n, x) RESULT(fn_val)

INTEGER, INTENT(IN)        :: n
REAL (dpv), INTENT(IN)      :: x
REAL (dpv)                  :: fn_val

fn_val = 0.5_dpv * LOG10(6.28_dpv*n) - n * LOG10(1.36_dpv*x/n)
RETURN
END FUNCTION envj

SUBROUTINE lpmns(m, n, x, pm, pd)

!    ========================================================
!    Purpose: Compute associated Legendre functions Pmn(x)
!             and Pmn'(x) for a given order
!    Input :  x --- Argument of Pmn(x)
!             m --- Order of Pmn(x),  m = 0,1,2,...,n
!             n --- Degree of Pmn(x), n = 0,1,2,...,N
!    Output:  PM(n) --- Pmn(x)
!             PD(n) --- Pmn'(x)
!    ========================================================

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dpv), INTENT(IN)   :: x
REAL (dpv), INTENT(OUT)  :: pm(0:n)
REAL (dpv), INTENT(OUT)  :: pd(0:n)

REAL (dpv)  :: pm0, pm1, pm2, pmk, x0
INTEGER    :: k

DO  k = 0, n
  pm(k) = 0.0_dpv
  pd(k) = 0.0_dpv
END DO
IF (ABS(x) == 1.0_dpv) THEN
  DO  k = 0, n
    IF (m == 0) THEN
      pm(k) = 1.0_dpv
      pd(k) = 0.5_dpv * k * (k+1)
      IF (x < 0.0) THEN
        pm(k) = (-1) ** k * pm(k)
        pd(k) = (-1) ** (k+1) * pd(k)
      END IF
    ELSE IF (m == 1) THEN
      pd(k) = 1.0D+300
    ELSE IF (m == 2) THEN
      pd(k) = -0.25_dpv * (k+2) * (k+1) * k * (k-1)
      IF (x < 0.0) pd(k) = (-1) ** (k+1) * pd(k)
    END IF
  END DO
  RETURN
END IF
x0 = ABS(1.0_dpv - x*x)
pm0 = 1.0_dpv
pmk = pm0
DO  k = 1, m
  pmk = (2*k-1) * SQRT(x0) * pm0
  pm0 = pmk
END DO
pm1 = (2*m+1) * x * pm0
pm(m) = pmk
pm(m+1) = pm1
DO  k = m + 2, n
  pm2 = ((2*k-1)*x*pm1 - (k+m-1)*pmk) / (k-m)
  pm(k) = pm2
  pmk = pm1
  pm1 = pm2
END DO
pd(0) = ((1-m)*pm(1) - x*pm(0)) / (x*x - 1.0)
DO  k = 1, n
  pd(k) = (k*x*pm(k)-(k+m)*pm(k-1)) / (x*x - 1.0_dpv)
END DO
RETURN
END SUBROUTINE lpmns
SUBROUTINE lpmn(mm, m, n, x, pm, pd)

!       =====================================================
!       Purpose: Compute the associated Legendre functions
!                Pmn(x) and their derivatives Pmn'(x)
!       Input :  x  --- Argument of Pmn(x)
!                m  --- Order of Pmn(x),  m = 0,1,2,...,n
!                n  --- Degree of Pmn(x), n = 0,1,2,...,N
!                mm --- Physical dimension of PM and PD
!       Output:  PM(m,n) --- Pmn(x)
!                PD(m,n) --- Pmn'(x)
!       =====================================================

INTEGER, INTENT(IN)     :: mm
INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
REAL (dpv), INTENT(IN)   :: x
REAL (dpv), INTENT(OUT)  :: pm(0:,0:)
REAL (dpv), INTENT(OUT)  :: pd(0:,0:)

REAL (dpv)  :: xq, xs
INTEGER    :: i, j, ls

DO  i = 0, n
  DO  j = 0, m
    pm(j,i) = 0.0D0
    pd(j,i) = 0.0D0
  END DO
END DO
pm(0,0) = 1.0D0
IF (ABS(x) == 1.0D0) THEN
  DO  i = 1, n
    pm(0,i) = x ** i
    pd(0,i) = 0.5D0 * i * (i+1.0D0) * x ** (i+1)
  END DO
  DO  j = 1, n
    DO  i = 1, m
      IF (i == 1) THEN
        pd(i,j) = 1.0D+300
      ELSE IF (i == 2) THEN
        pd(i,j) = -0.25D0 * (j+2) * (j+1) * j * (j-1) * x ** (j+1)
      END IF
    END DO
  END DO
  RETURN
END IF
ls = 1
IF (ABS(x) > 1.0D0) ls = -1
xq = SQRT(ls*(1.0D0 - x*x))
xs = ls * (1.0D0 - x*x)
DO  i = 1, m
  pm(i,i) = -ls * (2*i-1) * xq * pm(i-1,i-1)
END DO
DO  i = 0, m
  pm(i,i+1) = (2*i+1) * x * pm(i,i)
END DO
DO  i = 0, m
  DO  j = i + 2, n
    pm(i,j) = ((2*j-1)*x*pm(i,j-1) - (i+j-1)*pm(i,j-2)) / (j-i)
  END DO
END DO
pd(0,0) = 0.0D0
DO  j = 1, n
  pd(0,j) = ls * j * (pm(0,j-1) - x*pm(0,j)) / xs
END DO
DO  i = 1, m
  DO  j = i, n
    pd(i,j) = ls * i * x * pm(i,j) / xs + (j+i) * (j-i+1) / xq * pm(i-1,j)
  END DO
END DO
RETURN
END SUBROUTINE lpmn
!-------------------------------------------------------------------------------

SUBROUTINE radfun(alpha,n,l,r,radf)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,l
	REAL(DP), INTENT(IN) :: r
	REAL(DP), INTENT(OUT) :: radf
	REAL(DP) :: x, fact1, fact2, alpha

!	r = r * alpha
!	x =2.d0 * r / n
	x = 2.d0 * alpha * r / n
!	write(6,*)n,l,r,alpha
	call laguerre(n-l-1, 2 * l +1, x, radf)
!	write(6,*)"cp1"
	call factorial(n-l-1,fact1)
	call factorial(n+l,fact2)
	radf = radf * dexp(-x/2.d0) * (x)**l * 2.d0 / n / n * dsqrt(fact1/fact2) * dsqrt(alpha*alpha*alpha)

END SUBROUTINE radfun
!------------------------------------------------------------------------------


END MODULE specfun

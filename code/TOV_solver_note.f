        program TOV_solver

        parameter (Neos=5042, N=10)
        double precision Ec,Pc,dR,X,Y(3),DYDX(3),YOUT(3),
     &                   CSEVAL,R,M,YR,k2,lambda,M_old,
     &                   Ms(Neos),Rs(Neos),lambdas(Neos),
     &                   P0(Neos),E0(Neos),B0(Neos),C0(Neos),D0(Neos),
     &                   E1(Neos),P1(Neos),B1(Neos),C1(Neos),D1(Neos)

        common /EOS_PE/ P0, E0, B0, C0, D0
        common /EOS_EP/ E1, P1, B1, C1, D1
*---
* (1) Read in the EOS:

        open (UNIT=10, FILE='./EOS_simulated.dat')
        do 100 i = 1, Neos
                read (10,*) E0(i), P0(i)
                E1(i) = E0(i)
                P1(i) = P0(i)
100     continue
                close (10)

        call CSPLIN (Neos, P0, E0, B0, C0, D0)
        call CSPLIN (Neos, E1, P1, B1, C1, D1)

* (2) Solve extended TOV eqn for M, R, and lambda:

        dR = 0.000005d0

        do 200 i = Neos, 1, -N

           Ec = E0(i)
           Pc = P0(i) 

              X = dR
           Y(1) = Pc
           Y(2) = Ec * dR**3/3.d0
           Y(3) = 2.d0
           call DERIVS (X, Y, DYDX)
           do 210 k = 1, 300000
              call RK4 (Y, DYDX, 3, X, dR, YOUT)
              if (YOUT(1) .LT. 0.d0) goto 215
                 X = X + dR
              Y(1) = YOUT(1)
              Y(2) = YOUT(2)
              Y(3) = YOUT(3)
              call DERIVS (X, Y, DYDX)
 210       continue
           print *, '... Loop 210 !'
           ii = i + N
           goto 290
           stop '... Loop 210 !'
 215        R = 17.456d0 * X    ! in km
            M = 11.820d0 * Y(2) ! in Msun
           YR = Y(3)            ! dim-less

           call k2eval (R, M, YR, k2)
           lambda = 2.d0 * k2 * R**5 / 3.d0 / 6.67384d0 
           lambda = lambda / 1000.d0 ! in 10^36 g cm^2 s^2

           Ms(i) = M
           Rs(i) = R 
           lambdas(i) = lambda 
           if (M .LT. 0.5d0) THEN 
              ii = i
              goto 290
           end if

 200    continue

* (3) Save the result

 290    open (UNIT=20, FILE='./MRlda.dat')
        M_old = -1.d0
        do 300 j = ii, Neos, N
           if (Ms(j) .LT. M_old) THEN
              goto 310
           ELSE 
              WRITE (20,*) Ms(j), Rs(j), lambdas(j) 
              M_old = Ms(j)
           end if 
 300    continue
 310    close (20)
*---
        print *
        print *, '==================================================='
        print *, 'Neutron Star Properties                            '
        print *, '---------------------------------------------------'
        print *, 'MRlda.dat: M[Msun], R[km], lambda[10^36 g cm^2 s^2]'
        print *, '==================================================='
        print *
*---
        end
*=======================================================================
        SUBROUTINE DERIVS (X, Y, DYDX)

        parameter (Neos=5042)
        double precision X,Y(3),DYDX(3),CSEVAL,EOS,dPdE,
     &                   P0(Neos),E0(Neos),B0(Neos),C0(Neos),D0(Neos),
     &                   E1(Neos),P1(Neos),B1(Neos),C1(Neos),D1(Neos)

        common /EOS_PE/ P0, E0, B0, C0, D0
        common /EOS_EP/ E1, P1, B1, C1, D1
*---
         EOS = CSEVAL(Neos, Y(1), 0, P0, E0, B0, C0, D0)
        dPdE = CSEVAL(Neos,  EOS, 1, E1, P1, B1, C1, D1)

        F = (X - (EOS - Y(1)) * X**3) / (X - 2.d0 * Y(2))
        Q = (5.d0 * EOS + 9.d0 * Y(1) 
     &    + (EOS + Y(1)) / dPdE - 6.d0 / X**2) * X / (X - 2.d0 * Y(2))
     &    - 4.d0 * ((Y(2) + Y(1) * X**3) / X / (X - 2.d0 * Y(2)))**2
        
        DYDX(1) = - (EOS + Y(1)) * (Y(2) + Y(1) * X**3) 
     &          / X / (X - 2.d0 * Y(2))
        DYDX(2) = EOS * X**2 
        DYDX(3) = - (Y(3)**2 + F * Y(3) + Q * X**2) / X
*---
        return
        end
*=======================================================================
        SUBROUTINE RK4 (Y, DYDX, N, X, H, YOUT)

        parameter (NMAX=10)
        integer    N
        double precision H, HH, H6, X, XH, Y, DYDX, YOUT, YT, DYT, DYM
        dimension Y(N), DYDX(N), YOUT(N), YT(NMAX), DYT(NMAX), DYM(NMAX)
*---
        HH = H * 0.5
        H6 = H / 6.
        XH = X + HH
*---
        do 11 I = 1, N
           YT(I) = Y(I) + HH * DYDX(I)
 11     continue
*---
        call DERIVS (XH, YT, DYT)
        do 12 I = 1, N
           YT(I) = Y(I) + HH * DYT(I)
 12     continue
*---
        call DERIVS (XH, YT, DYM)
        do 13 I = 1, N
           YT(I)  = Y(I) + H * DYM(I)
           DYM(I) = DYT(I) + DYM(I)
 13     continue
*---
        call DERIVS (X+H, YT, DYT)
        do 14 I = 1, N
           YOUT(I) = Y(I) + H6 * (DYDX(I) + DYT(I) + 2. * DYM(I))
 14     continue
*---
        return
        end
*=======================================================================
        SUBROUTINE k2eval (R, M, YR, k2)

        double precision R,M,YR,k2,RsR,square,round,factor
*---
        RsR = 2.954d0 * M / R

        square = 26.d0 
     &         - 22.d0 * YR 
     &         + (3.d0 * YR - 2.d0) * RsR 
     &         + (1.d0 + YR) * RsR**2

        round = 6.d0
     &        - 3.d0 * YR
     &        + 1.5d0 * (5.d0 * YR - 8.d0) * RsR
     &        + 0.25d0 * square * RsR**2

        factor = (1.d0 - RsR)**2 * (2.d0 - YR + (YR - 1.d0) * RsR)

        k2 = 0.05d0 * RsR**5 * factor
     &     / (RsR * round + 3.d0 * factor * LOG(1.d0 - RsR))
*---
        return
        end
*=======================================================================
        SUBROUTINE CSPLIN (N, X, Y, B, C, D)
        integer N
        double precision X(N), Y(N), B(N), C(N), D(N)
C
C  THIS IS A MODifIED VERSION OF SPLINE, DESCRIBED IN THE NOTES
C  BY FORSYTHE MALCOLM AND MOLER
C
C  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
C  FOR A CUBIC INTERPOLATING SPLINE
C
C    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C    FOR  X(I) .LE. X .LE. X(I+1)
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
C    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
C    Y = THE ORDINATES OF THE KNOTS
C
C  OUTPUT..
C
C    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
C
C  USING  P  TO DENOTE DifFERENTIATION,
C
C    Y(I) = S(X(I))
C    B(I) = SP(X(I))
C    C(I) = SPP(X(I))/2
C    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
C
C  THE ACCOMPANYING FUNCTION SUBprogram  CSEVAL  CAN BE USED
C  TO EVALUATE THE SPLINE, ITS DERIVATIVE OR EVEN ITS 2ND DERIVATIVE.
C
C       SEE COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS BY
C           FORSYTHE, MALCOLM, AND MOLER FOR DETAILS
C
        integer LOUT, NM1, IB, I
        double precision T
        DATA LOUT/6/
C
C         CHECK INPUT FOR CONSISTENCY
C
        if (N .GE. 2) GO TO 1
        WRITE(LOUT, 1000)N
        return
C
    1   NM1 = N-1
        do 3 I = 1, NM1
           if (X(I) .LT. X(I+1)) GO TO 3
           print *, 'Wrong ordering at', I 
           WRITE(LOUT, 1001)
           return
    3   continue
C
    5   if ( N .EQ. 2 ) GO TO 50
C
C  SET UP TRIDIAGONAL SYSTEM
C
C  B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
C
        D(1) = X(2) - X(1)
        C(2) = (Y(2) - Y(1))/D(1)
        do 10 I = 2, NM1
           D(I) = X(I+1) - X(I)
           B(I) = 2.*(D(I-1) + D(I))
           C(I+1) = (Y(I+1) - Y(I))/D(I)
           C(I) = C(I+1) - C(I)
   10   continue
C
C  end CONDITIONS.  THIRD DERIVATIVES AT  X(1)  AND  X(N)
C  OBTAINED FROM DIVIDED DifFERENCES
C
        B(1) = -D(1)
        B(N) = -D(N-1)
        C(1) = 0.
        C(N) = 0.
        if ( N .EQ. 3 ) GO TO 15
        C(1) = C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1))
        C(N) = C(N-1)/(X(N)-X(N-2)) - C(N-2)/(X(N-1)-X(N-3))
        C(1) = C(1)*D(1)**2/(X(4)-X(1))
        C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
C
C  FORWARD ELIMINATION
C
   15   do 20 I = 2, N
           T = D(I-1)/B(I-1)
           B(I) = B(I) - T*D(I-1)
           C(I) = C(I) - T*C(I-1)
   20   continue
C
C  BACK SUBSTITUTION
C
        C(N) = C(N)/B(N)
        do 30 IB = 1, NM1
           I = N-IB
           C(I) = (C(I) - D(I)*C(I+1))/B(I)
   30   continue
C
C  C(I) IS NOW THE SIGMA(I) OF THE TEXT
C
C  COMPUTE POLYNOMIAL COEFFICIENTS
C
        B(N) = (Y(N) - Y(NM1))/D(NM1) + D(NM1)*(C(NM1) + 2.*C(N))
        do 40 I = 1, NM1
           B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
           D(I) = (C(I+1) - C(I))/D(I)
           C(I) = 3.*C(I)
   40   continue
        C(N) = 3.*C(N)
        D(N) = D(N-1)
        return
   50   B(1) = (Y(2)-Y(1))/(X(2)-X(1))
        C(1) = 0.
        D(1) = 0.
        return
 1000   format('-N < 2 IN CSPLIN call--',I10)
 1001   format('-X IS NOT IN ASCendING ORDER IN CSPLIN call')
        end
        function CSEVAL (N, U, IDERIV, X, Y, B, C, D)
        integer N, IDERIV
        double precision U, X(N), Y(N), B(N), C(N), D(N), CSEVAL
C
C  THIS IS A MODifIED VERSION OF SEVAL, DESCRIBED IN THE NOTES
C  BY FORSYTHE MALCOLM AND MOLER
C
C  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION, ITS FIRST
C  DERIVATIVE OR SECOND DERIVATIVE, NAMELY:
C
C   IDERIV=0:    Y(I)+B(I)*(U-X(I))+C(I)*(U-X(I))**2+D(I)*(U-X(I))**3
C   IDERIV=1:    B(I)+2*C(I)*(U-X(I))+3*D(I)*(U-X(I))**2
C   IDERIV=2:    2*C(I)+6*D(I)*(U-X(I))
C
C    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
C
C  if  U .LT. X(1) THEN  I = 1  IS USED.
C  if  U .GE. X(N) THEN  I = N  IS USED.
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS
C    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
C    IDERIV = 0 TO EVALUATE S(U)
C           = 1 TO EVALUATE SP(U)
C           = 2 TO EVALUATE SPP(U)
C    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
C    B,C,D = ARRAYS OF SPLINES COEFFICIENTS COMPUTED BY SUBROUTINE CSPLI
C
C  if  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS call, THEN A
C    BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
C
        integer I, J, K, LOUT
        double precision DX
        DATA I, LOUT/1, 6/
C
        if (IDERIV .GE. 0 .AND. IDERIV .LE. 2) GO TO 5
        WRITE(LOUT, 1000)IDERIV
        return
C
    5   if ( I .GE. N ) I = 1
        if ( U .LT. X(I) ) GO TO 10
        if ( U .LE. X(I+1) ) GO TO 30
C
C  BINARY SEARCH
C
   10   I = 1
        J = N + 1
   20   K = (I+J)/2
        if ( U .LT. X(K) ) J = K
        if ( U .GE. X(K) ) I = K
        if ( J .GT. I+1 ) GO TO 20
C
C  EVALUATE SPLINE
C
   30   DX = U - X(I)
        J = IDERIV+1
        GO TO (40,50,60), J
C  COMPUTE S(X)
   40   CSEVAL = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
        return
C  COMPUTE SP(X)
   50   CSEVAL = B(I) + DX*(2*C(I) + 3*DX*D(I))
        return
C  COMPUTE SPP(X)
   60   CSEVAL = 2*C(I) + 6*DX*D(I)
        return
 1000   FORMAT('-IDERIV IS INVALID IN CSEVAL call--', I10)
        end
C=======================================================================


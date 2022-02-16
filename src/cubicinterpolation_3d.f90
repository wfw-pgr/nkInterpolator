! ====================================================== !
! === 3Dcubic interpolation ( Lagrange )             === !
! ====================================================== !
subroutine cubicinterpolation_3d( xRef, xItp, LI, LJ, LK, nItp )
  implicit none

  integer         , intent(in)    :: LI, LJ, LK, nItp
  double precision, intent(in)    :: xRef(4,LI,LJ,LK)
  double precision, intent(inout) :: xItp(4,nItp)
  integer         , parameter     :: lun = 50
  integer         , parameter     :: x1_=1, x2_=2, x3_=3, x0_=4
  integer                         :: i, j, k, m, ki, ip, jp, kp
  integer         , allocatable   :: xItp_index(:,:)
  double precision, allocatable   :: BMatrix(:,:), DMatrix(:,:), AMatrix(:,:)
  double precision                :: x1Min, x2Min, x3Min, dx1Inv, dx2Inv, dx3Inv, fi
  double precision, allocatable   :: xItp_norm(:), avec(:,:,:), pvec(:,:,:)
  logical                         :: Flag__OutOfIndex
  logical, parameter              :: Flag__WarningDisplay    = .true.
  logical, parameter              :: Flag__saveMatrices      = .false.
  logical, parameter              :: Flag__StoreZeroValue    = .false.
  logical, parameter              :: Flag__StoreNearestPoint = .true.
  logical, parameter              :: Flag__LinearInterpolate = .false.


  ! ------------------------------------------------------ !
  ! --- [1] Prepare Bmatrix Coefficients               --- !
  ! ------------------------------------------------------ !
  allocate( BMatrix(64,64), DMatrix(64,64), AMatrix(64,64) )
  call generate__Matrices( BMatrix, DMatrix, AMatrix )
  ! ------------------------------------------------------ !
  ! --- [2] interpolation pt. ==> grid number          --- !
  ! ------------------------------------------------------ !
  allocate( xItp_index(3,nItp) )
  x1Min  = xRef(x1_,1,1,1)
  x2Min  = xRef(x2_,1,1,1)
  x3Min  = xRef(x3_,1,1,1)
  dx1Inv = 1.d0 / ( xRef(x1_,2,1,1) - xRef(x1_,1,1,1) )
  dx2Inv = 1.d0 / ( xRef(x2_,1,2,1) - xRef(x2_,1,1,1) )
  dx3Inv = 1.d0 / ( xRef(x3_,1,1,2) - xRef(x3_,1,1,1) )
  do ki=1, nItp
     xItp_index(x1_,ki) = ( ceiling( ( xItp(x1_,ki)-x1Min )*dx1Inv ) - 1 ) + 1
     xItp_index(x2_,ki) = ( ceiling( ( xItp(x2_,ki)-x2Min )*dx2Inv ) - 1 ) + 1
     xItp_index(x3_,ki) = ( ceiling( ( xItp(x3_,ki)-x3Min )*dx3Inv ) - 1 ) + 1
  enddo

  ! ------------------------------------------------------ !
  ! --- [3] cubic interpolation loop                   --- !
  ! ------------------------------------------------------ !
  allocate( xItp_norm(3), avec(4,4,4), pvec(4,4,4) )
  do ki=1, nItp
     !  -- [3-1] out-of-index exception                --  !
     Flag__OutOfIndex = .false.
     if ( ( xItp_index(x1_,ki).le.1 ).or.( xItp_index(x1_,ki).ge.LI-1 ) ) Flag__OutOfIndex = .true.
     if ( ( xItp_index(x2_,ki).le.1 ).or.( xItp_index(x2_,ki).ge.LJ-1 ) ) Flag__OutOfIndex = .true.
     if ( ( xItp_index(x3_,ki).le.1 ).or.( xItp_index(x3_,ki).ge.LK-1 ) ) Flag__OutOfIndex = .true.
     if ( Flag__OutOfIndex ) then
        if ( Flag__StoreNearestPoint ) then
           ip = max( min( xItp_index(x1_,ki), LI ), 1 )
           jp = max( min( xItp_index(x2_,ki), LJ ), 1 )
           kp = max( min( xItp_index(x3_,ki), LK ), 1 )
           xItp(x0_,ki) = xRef( x0_,ip,jp,kp )
           if ( Flag__WarningDisplay ) then
              write(6,*) '[cubInterpMod] WARNING!! out-of-index Exception!!  :: (ip,jp,kp) == ', ip, jp, kp
           endif
           cycle
        endif
        if ( Flag__StoreZeroValue ) then
           xItp(x0_,ki) = 0.d0
           if ( Flag__WarningDisplay ) then
              write(6,*) '[cubInterpMod] WARNING!! out-of-index Exception!!  :: (ip,jp,kp) == ', ip, jp, kp
           endif
           cycle
        endif
        ! if ( Flag__LinearInterpolate ) then
        !    if ( xItp_index(x1_,ki).le.1 ) then
        !    ! xItp(x0_,ki) = 
        !    if ( Flag__WarningDisplay ) then
        !       write(6,*) '[cubInterpMod] WARNING!! out-of-index Exception!!  :: (ip,jp,kp) == ', ip, jp, kp
        !    endif
        !    cycle
        ! endif
     endif
     do k=1, 4
        do j=1, 4
           do i=1, 4
              avec(i,j,k) = 0.d0
              pvec(i,j,k) = xRef( x0_, xItp_index(x1_,ki)-2+i, xItp_index(x2_,ki)-2+j, xItp_index(x3_,ki)-2+k )
           enddo
        enddo
     enddo
     !  -- [3-2] Coefficient calculation               --  !
     call MatrixVectorMultiply( AMatrix, pvec, avec, 64, 64 )
     !  -- [3-3] polynomial interpolation              --  !
     ip             = xItp_index(x1_,ki)
     jp             = xItp_index(x2_,ki)
     kp             = xItp_index(x3_,ki)
     xItp_norm(x1_) = ( xItp(x1_,ki) - xRef(x1_,ip,jp,kp) )*dx1Inv
     xItp_norm(x2_) = ( xItp(x2_,ki) - xRef(x2_,ip,jp,kp) )*dx2Inv
     xItp_norm(x3_) = ( xItp(x3_,ki) - xRef(x3_,ip,jp,kp) )*dx3Inv
     fi             = 0.d0
     do k=1, 4
        do j=1, 4
           do i=1, 4
              fi = fi + avec(i,j,k)*xItp_norm(x1_)**(i-1) * xItp_norm(x2_)**(j-1) * xItp_norm(x3_)**(k-1)
           enddo
        enddo
     enddo
     xItp(x0_,ki) = fi
  enddo
  ! ------------------------------------------------------ !
  ! --- [4] post process                               --- !
  ! ------------------------------------------------------ !
  deallocate( BMatrix, DMatrix, AMatrix )
  deallocate( xItp_index, xItp_norm, avec, pvec )
  return
  
contains

  ! ====================================================== !
  ! === generate ( B, D, A ) Matrices                  === !
  ! ====================================================== !
  subroutine generate__Matrices( BMatrix, DMatrix, AMatrix )
    implicit none
    double precision, intent(out) :: BMatrix(64,64), DMatrix(64,64), AMatrix(64,64)

    call generate__BMatrix( BMatrix )
    call generate__DMatrix( DMatrix )
    call generate__AMatrix( BMatrix, DMatrix, AMatrix )
    if ( Flag__saveMatrices ) then
       call save__Matrices( BMatrix, DMatrix, AMatrix )
    endif
    
    return
  end subroutine generate__Matrices
  

  ! ====================================================== !
  ! === generate Interpolation (B) Matrix              === !
  ! ====================================================== !
  subroutine generate__BMatrix( BMatrix )
    implicit none
    integer                       :: i, j, k, idx, ex, ey, ez
    integer                       :: exyz(64,3)
    double precision              :: xn, yn, zn
    double precision              :: corners(8,3)
    double precision, intent(out) :: BMatrix(64,64)

    ! ------------------------------------------------------ !
    ! --- [1] index preparation                          --- !
    ! ------------------------------------------------------ !
    !  -- [1-1] normalized corner points [ 0.0, 1.0 ]    --  !
    do k=1, 2
       do j=1, 2
          do i=1, 2
             idx            = 4*(k-1) + 2*(j-1) + (i-1) + 1
             corners(idx,1) = i-1
             corners(idx,2) = j-1
             corners(idx,3) = k-1
          enddo
       enddo
    enddo
    !  -- [1-2] exponential number x^e                   --  !
    do k=1, 4
       do j=1, 4
          do i=1, 4
             idx         = 16*(k-1) + 4*(j-1) + (i-1) + 1
             exyz(idx,1) = i-1
             exyz(idx,2) = j-1
             exyz(idx,3) = k-1
          enddo
       enddo
    enddo
    ! ------------------------------------------------------ !
    ! --- [2] BMatrix Calculation                        --- !
    ! ------------------------------------------------------ !
    do i=1, 64
       ex = exyz(i,1)
       ey = exyz(i,2)
       ez = exyz(i,3)
       do k=1, 8
          xn               = corners(k,1)
          yn               = corners(k,2)
          zn               = corners(k,3)
          BMatrix(0*8+k,i) = xn**ex             * yn**ey             * zn**ez
          BMatrix(1*8+k,i) = ex*xn**(abs(ex-1)) * yn**ey             * zn**ez
          BMatrix(2*8+k,i) = xn**ex             * ey*yn**(abs(ey-1)) * zn**ez
          BMatrix(3*8+k,i) = xn**ex             * yn**ey             * ez*zn**(abs(ez-1))
          BMatrix(4*8+k,i) = ex*xn**(abs(ex-1)) * ey*yn**(abs(ey-1)) * zn**ez
          BMatrix(5*8+k,i) = ex*xn**(abs(ex-1)) * yn**ey             * ez*zn**(abs(ez-1))
          BMatrix(6*8+k,i) = xn**ex             * ey*yn**(abs(ey-1)) * ez*zn**(abs(ez-1))
          BMatrix(7*8+k,i) = ex*xn**(abs(ex-1)) * ey*yn**(abs(ey-1)) * ez*zn**(abs(ez-1))
       enddo
    enddo
    return
  end subroutine generate__BMatrix


  ! ====================================================== !
  ! ===  generate Differential Matrix ( D )            === !
  ! ====================================================== !
  subroutine generate__DMatrix( DMatrix )
    implicit none
    integer                       :: i, j, m
    integer                       :: cuboid(8)
    double precision, intent(out) :: DMatrix(64,64)

    ! ------------------------------------------------------ !
    ! --- [1] preparation cboid / DMatrix                --- !
    ! ------------------------------------------------------ !
    data cuboid/22,23,26,27,38,39,42,43/
    do j=1, 64
       do i=1, 64
          DMatrix(i,j) = 0.d0
       enddo
    enddo
    
    ! ------------------------------------------------------ !
    ! --- [1] preparation 'cuboid' Coordinates           --- !
    ! ------------------------------------------------------ !
    do m=1, 8
       DMatrix( m+0 , cuboid(m)    ) = 1.d0
    enddo
    do m=1, 8
       DMatrix( m+8 , cuboid(m)-1  ) = -0.5d0
       DMatrix( m+8 , cuboid(m)+1  ) = +0.5d0
    enddo
    do m=1, 8
       DMatrix( m+16, cuboid(m)-4  ) = -0.5d0
       DMatrix( m+16, cuboid(m)+4  ) = +0.5d0
    enddo
    do m=1, 8
       DMatrix( m+24, cuboid(m)-16 ) = -0.5d0
       DMatrix( m+24, cuboid(m)+16 ) = +0.5d0
    enddo
    do m=1, 8
       DMatrix( m+32, cuboid(m)+5  ) = +0.25d0
       DMatrix( m+32, cuboid(m)-3  ) = -0.25d0
       DMatrix( m+32, cuboid(m)+3  ) = -0.25d0
       DMatrix( m+32, cuboid(m)-5  ) = +0.25d0
    enddo
    do m=1, 8
       DMatrix( m+40, cuboid(m)+17 ) = +0.25d0
       DMatrix( m+40, cuboid(m)-15 ) = -0.25d0
       DMatrix( m+40, cuboid(m)+15 ) = -0.25d0
       DMatrix( m+40, cuboid(m)-17 ) = +0.25d0
    enddo
    do m=1, 8
       DMatrix( m+48, cuboid(m)+20 ) = +0.25d0
       DMatrix( m+48, cuboid(m)-12 ) = -0.25d0
       DMatrix( m+48, cuboid(m)+12 ) = -0.25d0
       DMatrix( m+48, cuboid(m)-20 ) = +0.25d0
    enddo
    do m=1, 8
       DMatrix( m+56, cuboid(m)+21 ) = +0.125d0
       DMatrix( m+56, cuboid(m)+13 ) = -0.125d0
       DMatrix( m+56, cuboid(m)+19 ) = -0.125d0
       DMatrix( m+56, cuboid(m)+11 ) = +0.125d0
       DMatrix( m+56, cuboid(m)-11 ) = -0.125d0
       DMatrix( m+56, cuboid(m)-19 ) = +0.125d0
       DMatrix( m+56, cuboid(m)-13 ) = +0.125d0
       DMatrix( m+56, cuboid(m)-21 ) = -0.125d0
    enddo
    return
  end subroutine generate__DMatrix


  ! ====================================================== !
  ! === generate united Matrix ( A )                   === !
  ! ====================================================== !
  subroutine generate__AMatrix( BMatrix, DMatrix, AMatrix )
    implicit none
    double precision, intent(in)  :: BMatrix(64,64), DMatrix(64,64)
    double precision, intent(out) :: AMatrix(64,64)
    double precision, allocatable :: CMatrix(:,:), work(:), ipiv(:)
    double precision              :: alpha, beta
    integer                       :: i, j, M, N, K, LDC, LDD, LDA, lwork, info

    ! ------------------------------------------------------ !
    ! --- [1] prepare Matrix                             --- !
    ! ------------------------------------------------------ !
    allocate( CMatrix(64,64) )
    do j=1, 64
       do i=1, 64
          CMatrix(i,j) = BMatrix(i,j)
          AMatrix(i,j) = 0.d0
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] Matrix Inversion                           --- !
    ! ------------------------------------------------------ !
    N      = 64
    LDC    = 64
    lwork  = N**2
    allocate( work(lwork), ipiv(N) )
    call dgetrf( N, N, CMatrix, LDC, ipiv,              info )
    call dgetri( N,    CMatrix, LDC, ipiv, work, lwork, info )

    ! ------------------------------------------------------ !
    ! --- [3] Matrix Multiplication                      --- !
    ! ------------------------------------------------------ !
    M      = 64
    N      = 64
    K      = 64
    LDC    = 64
    LDD    = 64
    LDA    = 64
    alpha  = 1.d0
    beta   = 0.d0
    call dgemm( 'N', 'N', M, N, K, alpha, CMatrix, LDC, DMatrix, LDD, beta, AMatrix, LDA )
    
    return
  end subroutine generate__AMatrix
    
  
  ! ====================================================== !
  ! === save__Matrices                                 === !
  ! ====================================================== !
  subroutine save__Matrices( BMatrix, DMatrix, AMatrix )
    implicit none
    integer         , parameter  :: lun = 50
    double precision, intent(in) :: BMatrix(64,64), DMatrix(64,64), AMatrix(64,64)
    integer                      :: i, j

    ! ------------------------------------------------------ !
    ! --- [1] save BMatrix                               --- !
    ! ------------------------------------------------------ !
    open(lun,file='BMatrix.dat',form='formatted')
    do i=1, 64
       write(lun,*) ( BMatrix(i,j), j=1, 64 )
    enddo
    close(lun)
    ! ------------------------------------------------------ !
    ! --- [2] save DMatrix                               --- !
    ! ------------------------------------------------------ !
    open(lun,file='DMatrix.dat',form='formatted')
    do i=1, 64
       write(lun,*) ( DMatrix(i,j), j=1, 64 )
    enddo
    close(lun)
    ! ------------------------------------------------------ !
    ! --- [3] save AMatrix                               --- !
    ! ------------------------------------------------------ !
    open(lun,file='AMatrix.dat',form='formatted')
    do i=1, 64
       write(lun,*) ( AMatrix(i,j), j=1, 64 )
    enddo
    close(lun)

    return
  end subroutine save__Matrices


  
  ! ====================================================== !
  ! === BLAS Matrix Vector Multiply Wrapper            === !
  ! ====================================================== !
  subroutine MatrixVectorMultiply( Amat, xvec, yvec, LM, LN )
    implicit none
    integer         , intent(in)  :: LM, LN
    double precision, intent(in)  :: Amat(LM,LN), xvec(LN)
    double precision, intent(out) :: yvec(LM)
    integer         , parameter   :: incX=1, incY=1
    double precision, parameter   :: c1=1.d0, c2=0.d0

    ! ------------------------------------------------------ !
    ! --- [1] call DGEMV                                 --- !
    ! ------------------------------------------------------ !
    yvec(:) = 0.d0
    call DGEMV( 'N', LM, LN, c1, Amat, LM, xvec, incX, c2, yvec, incY )
    
    return
  end subroutine MatrixVectorMultiply


  ! ! ====================================================== !
  ! ! === linear interpolation at edge                   === !
  ! ! ====================================================== !
  ! subroutine interpolate__linear
  !   implicit none
  !   integer         , intent(in)  :: 
  !   double precision, intent(out) :: 
  !   if ( i_index.le.1 ) then
  !      ( xItp(x_,ip) - xRef(x_,1,jk,kk) ) * ( xRef(v_,2,jk,kk) - xRef(v_,1,jk,kk) ) / ( xRef(x_,2,jk,kk) - xRef(x_,1,jk,kk) )
  !      xRef(x_,)

  !      .or.( i_index.ge.LI-1 ) ) Flag__OutOfIndex = .true.
  !   if ( ( xItp_index(x2_,ki).le.1 ).or.( xItp_index(x2_,ki).ge.LJ-1 ) ) Flag__OutOfIndex = .true.
  !   if ( ( xItp_index(x3_,ki).le.1 ).or.( xItp_index(x3_,ki).ge.LK-1 ) ) Flag__OutOfIndex = .true.
    

  !   return
  ! end subroutine interpolate__linear

  
    

end subroutine cubicinterpolation_3d

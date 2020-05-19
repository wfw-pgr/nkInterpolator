! ====================================================== !
! === handle all of the points to be interpolated    === !
! ====================================================== !
subroutine barycentric__interpolator( nodes , points , simplex, pwhere, &
     &                                nNodes, nPoints, nSimplex )
  implicit none
  integer         , intent(in)    :: nNodes, nSimplex, nPoints
  double precision, intent(in)    :: nodes  (3,nNodes  )
  double precision, intent(inout) :: points (3,nPoints )
  integer         , intent(in)    :: simplex(3,nSimplex), pwhere(nPoints)
  integer                         :: ipt, i1, i2, i3
  
  do ipt=1, nPoints
     
     if ( pwhere(ipt).eq.-1 ) then
        call find__nearestpoint( nodes, points(:,ipt), nNodes )
     else
        ! -- simplex, pwhere :: begin from 0 -- !
        i1 = simplex(1,pwhere(ipt)+1 ) + 1
        i2 = simplex(2,pwhere(ipt)+1 ) + 1
        i3 = simplex(3,pwhere(ipt)+1 ) + 1
        call barycentric__interpolate( nodes(:,i1), nodes(:,i2), nodes(:,i3), points(:,ipt) )
     endif
     
  enddo
  
  return
contains

  ! ====================================================== !
  ! === interpolation in barycentric coordinate        === !
  ! ====================================================== !
  subroutine barycentric__interpolate( v1, v2, v3, xp )
    implicit none
    double precision, intent(in)    :: v1(3), v2(3), v3(3)
    double precision, intent(inout) :: xp(3)
    double precision                :: denom, pI, w1, w2, w3
    integer         , parameter     :: x_=1, y_=2, z_=3

    denom = ( ( v2(y_)-v3(y_) )*( v1(x_)-v3(x_) ) + ( v3(x_)-v2(x_) )*( v1(y_)-v3(y_) ) )
    if ( denom.eq.0.d0 ) then
       pI = 0.d0
    else
       pI = 1.d0 / denom
    endif
    w1    = pI * ( ( v2(y_)-v3(y_) )*( xp(x_)-v3(x_) ) + ( v3(x_)-v2(x_) )*( xp(y_)-v3(y_) ) )
    w2    = pI * ( ( v3(y_)-v1(y_) )*( xp(x_)-v3(x_) ) + ( v1(x_)-v3(x_) )*( xp(y_)-v3(y_) ) )
    w3    = 1.d0  - w1 - w2

    xp(z_) = w1*v1(z_) + w2*v2(z_) + w3*v3(z_)

    return
  end subroutine barycentric__interpolate


  ! ====================================================== !
  ! === find nearest node                              === !
  ! ====================================================== !
  subroutine find__nearestpoint( nodes, point, nNodes )
    implicit none
    integer         , intent(in)    :: nNodes
    double precision, intent(in)    :: nodes(3,nNodes)
    double precision, intent(inout) :: point(3)
    integer         , parameter     :: x_=1, y_=2, z_=3
    integer                         :: iN  , i_min
    double precision                :: dist, d_min

    ! --  initial setup  -- !
    iN    = 1
    i_min = iN
    d_min = sqrt( ( point(x_) - nodes(x_,iN) )**2 + ( point(y_) - nodes(x_,iN) )**2 )
    ! --  find nearest   -- !
    do iN=2, nNodes
       dist = sqrt( ( point(x_) - nodes(x_,iN) )**2 + ( point(y_) - nodes(x_,iN) )**2 )
       if ( dist.lt.d_min ) then
          d_min = dist
          i_min = iN
       endif
    enddo
    ! --  substitute     -- !
    point(z_) = nodes(z_,i_min)

    return
  end subroutine find__nearestpoint
  
end subroutine barycentric__interpolator

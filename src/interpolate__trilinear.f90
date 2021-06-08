

! ====================================================== !
! === trilinear interpolation                        === !
! ====================================================== !
subroutine interpolate__trilinear( gData, pData, xMin, yMin, zMin, dx, dy, dz, &
     & LI, LJ, LK, nData, i_outOfRangeMode )
  implicit none
  integer         , intent(in)    :: LI, LJ, LK, nData, i_outOfRangeMode
  double precision, intent(in)    :: dx, dy, dz, xMin, yMin, zMin
  double precision, intent(in)    :: gData(LI,LJ,LK)
  double precision, intent(inout) :: pData(4,nData)
  integer                         :: ip, jp, kp, ii, jj, kk, ir, jr, kr, iData
  double precision                :: xp, yp, zp, dxInv, dyInv, dzInv, sum, xMax, yMax, zMax
  double precision                :: sfx(-2:2), sfy(-2:2), sfz(-2:2)
  integer         , parameter     :: x_=1, y_=2, z_=3, v_=4
  double precision                :: fill_value
  double precision, parameter     :: AlertLargeValue = 1.d10
  logical                         :: Flag__FillValue = .false.
  logical                         :: Flag__MinMax    = .false.

  ! ------------------------------------------------------ !
  ! --- [1] preparation                                --- !
  ! ------------------------------------------------------ !
  xMax  = dx*dble(LI-1) + xMin
  yMax  = dy*dble(LJ-1) + yMin
  zMax  = dz*dble(LK-1) + zMin
  dxInv = 1.d0 / dx
  dyInv = 1.d0 / dy
  dzInv = 1.d0 / dz

  if      ( i_outOfRangeMode.eq.0 ) then
     Flag__FillValue = .true.
     fill_value      = 0.d0
  else if ( i_outOfRangeMode.eq.1 ) then
     Flag__MinMax    = .true.
  else if ( i_outOfRangeMode.eq.2 ) then
     Flag__FillValue = .true.
     fill_value      = AlertLargeValue
  else
     write(6,*) "[interpolate__trilinear.f90] illegal i_outOfRangeMode !! :: ", i_outOfRangeMode
     stop
  endif
  
  ! ------------------------------------------------------ !
  ! --- [2] forced inside Grid Data mode               --- !
  ! ------------------------------------------------------ !
  if ( Flag__MinMax ) then
     do iData=1, nData
        xp   = min( max( pData(x_,iData)-xMin, 0.d0 ), xMax-xMin )
        yp   = min( max( pData(y_,iData)-yMin, 0.d0 ), yMax-yMin )
        zp   = min( max( pData(z_,iData)-zMin, 0.d0 ), zMax-zMin )
        ip   = nint( xp*dxInv )
        jp   = nint( yp*dyInv )
        kp   = nint( zp*dzInv )
        sfx  = shapef( xp, ip, dxInv )
        sfy  = shapef( yp, jp, dyInv )
        sfz  = shapef( zp, kp, dzInv )
        sum  = 0.d0
        do kk=-2, 2
           do jj=-2, 2
              do ii=-2, 2
                 ir  = ip+1 + ii
                 jr  = jp+1 + jj
                 kr  = kp+1 + kk
                 ir  = max( min( ir, LI ), 1 )
                 jr  = max( min( jr, LJ ), 1 )
                 kr  = max( min( kr, LK ), 1 )
                 sum = sum + sfx(ii)*sfy(jj)*sfz(kk) * gData(ir,jr,kr)
              enddo
           enddo
        enddo
        pData(v_,iData) = sum
     enddo
  end if

  ! ------------------------------------------------------ !
  ! --- [3] Fill value for out of region ver.          --- !
  ! ------------------------------------------------------ !
  if ( Flag__FillValue ) then
     
     do iData=1, nData

        if (   ( pData(x_,iData).lt.xMin ).or.( pData(x_,iData).gt.xMax ).or. &
             & ( pData(y_,iData).lt.yMin ).or.( pData(y_,iData).gt.yMax ).or. &
             & ( pData(z_,iData).lt.zMin ).or.( pData(z_,iData).gt.zMax ) ) then
           pData(v_,iData) = fill_value
        else
           
           xp   = pData(x_,iData)-xMin
           yp   = pData(y_,iData)-yMin
           zp   = pData(z_,iData)-zMin
           ip   = nint( xp*dxInv )
           jp   = nint( yp*dyInv )
           kp   = nint( zp*dzInv )
           sfx  = shapef( xp, ip, dxInv )
           sfy  = shapef( yp, jp, dyInv )
           sfz  = shapef( zp, kp, dzInv )
           sum  = 0.d0
           do kk=-2, 2
              do jj=-2, 2
                 do ii=-2, 2
                    ir  = ip+1 + ii
                    jr  = jp+1 + jj
                    kr  = kp+1 + kk
                    ir  = max( min( ir, LI ), 1 )
                    jr  = max( min( jr, LJ ), 1 )
                    kr  = max( min( kr, LK ), 1 )
                    sum = sum + sfx(ii)*sfy(jj)*sfz(kk) * gData(ir,jr,kr)
                 enddo
              enddo
           enddo
           pData(v_,iData) = sum

        endif
     enddo
  end if

  return
contains
  
  ! ====================================================== !
  ! === shape Function                                 === !
  ! ====================================================== !
  function shapef( rp, lp, drInv )
    implicit none
    integer         , intent(in) :: lp
    double precision, intent(in) :: rp, drInv
    double precision             :: delta
    double precision             :: shapef(-2:+2)

    delta     = rp*drInv - dble(lp)
    shapef(-2) = 0.d0
    shapef(-1) = 0.d0
    shapef( 0) = 1.d0 - abs( delta )
    shapef(+1) = 0.d0
    shapef(+2) = 0.d0
    shapef( int( sign( 1.d0, delta ) ) ) = abs( delta )

    return
  end function shapef

end subroutine interpolate__trilinear

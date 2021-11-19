import math, os, sys
import numpy         as np

# ========================================================= #
# ===  interpolate__polynomial                          === #
# ========================================================= #

def interpolate__polynomial( xp=None, xd=None, fx=None, fpx=None, fppx=None, \
                             order=5, nDerivative=1 ):

    # ------------------------------------------------- #
    # --- [1] size check                            --- #
    # ------------------------------------------------- #
    nCoef   = order       + 1
    nEqType = nDerivative + 1
    nPoint  = xd.shape[0]

    if ( nPoint * nEqType != nCoef ):
        print( "[interpolate__polynomial.py] information is not sufficient.... [ERROR]" )
        print( "[interpolate__polynomial.py] order       ::  {0}".format( order       ) )
        print( "[interpolate__polynomial.py]  => nCoef   ::  {0}".format( nCoef       ) )
        print()
        print( "[interpolate__polynomial.py] nPoint      ::  {0}".format( nPoint      ) )
        print()
        print( "[interpolate__polynomial.py] nDerivative ::  {0}".format( nDerivative ) )
        print( "[interpolate__polynomial.py]  => nEqType ::  {0}".format( nEqType     ) )
        print( "[interpolate__polynomial.py] nCoef != nEqType * nPoint !!! [ERROR]" )
        sys.exit()
        
    # ------------------------------------------------- #
    # --- [2] matrix making                         --- #
    # ------------------------------------------------- #
    Amat    =  np.zeros( (nCoef,nCoef) )
    if   ( nDerivative == 0 ):
        lhs = np.copy( fx )
    elif ( nDerivative == 1 ):
        lhs = np.concatenate( [fx,fpx] )
    elif ( nDerivative == 2 ):
        lhs = np.concatenate( [fx,fpx,fppx] )

    if ( lhs.shape[0] != nCoef ):
        print( "[interpolate__polynomial.py] data        is not sufficient.... [ERROR]" )
        if ( nDerivative >= 0 ):
            print( "[interpolate__polynomial.py] xd.shape   ::  {0}".format( xd.shape   ) )
            print( "[interpolate__polynomial.py] fx.shape   ::  {0}".format( fx.shape   ) )
        if ( nDerivative >= 1 ):
            print( "[interpolate__polynomial.py] fpx.shape  ::  {0}".format( fpx.shape  ) )
        if ( nDerivative >= 2 ):
            print( "[interpolate__polynomial.py] fppx.shape ::  {0}".format( fppx.shape ) )
        

    for iDeriv in range( nEqType ):
        irow  = nPoint * iDeriv
        for ik in range( 0, nPoint ):
            for jk in range( iDeriv, nCoef ):
                deriv_coef       = math.factorial( jk ) / math.factorial( jk-iDeriv )
                Amat[irow+ik,jk] = deriv_coef * xd[ik]**( float(jk-iDeriv) )
            
    # ------------------------------------------------- #
    # --- [3] solve coefficient                     --- #
    # ------------------------------------------------- #
    coef   = np.dot( np.linalg.inv( Amat ), lhs )

    # ------------------------------------------------- #
    # --- [4] interpolation                         --- #
    # ------------------------------------------------- #
    ret    = np.zeros( (xp.shape[0],) )
    for ik in range( nCoef ):
        ret[:] = ret[:] + coef[ik]*xp[:]**( float(ik) )
    return( ret )


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):

    # ------------------------------------------------- #
    # --- [1] test profile generation               --- #
    # ------------------------------------------------- #
    def arbitral__func( xin, coef=None, order=4 ):
        if ( coef is None ):
            coef = np.random.uniform( size=order+1 )
            print( coef )
        ret = np.zeros( (xin.shape[0],) )
        for ik in range( order+1 ):
            ret[:] = ret[:] + coef[ik]*xin[:]**( float(ik) )
        return( ret )

    xp = np.linspace( 0.0, 2.4, 101 )
    yp = arbitral__func( xp )


    # ------------------------------------------------- #
    # --- [2] case1: use ~1st derivatives and 3 pts --- #
    # ------------------------------------------------- #
    # ipt = xp.shape[0] // 2
    # xd   = np.array( [ xp[1], xp[ipt], xp[-2] ] )
    # fx   = np.array( [ yp[1], yp[ipt], yp[-2] ] )
    # fpx  = np.array( [ ( yp[    2] - yp[    0]) / ( xp[    2] - xp[    0] ), \
    #                    ( yp[ipt+1] - yp[ipt-1]) / ( xp[ipt+1] - xp[ipt-1] ), \
    #                    ( yp[   -1] - yp[   -3]) / ( xp[   -1] - xp[   -3] ) ] )
    # fppx = np.array( [ ( yp[    2] - 2*yp[  1] + yp[    0] ) / ( 0.5*( xp[    2] - xp[    0] ) )**2, \
    #                    ( yp[ipt+1] - 2*yp[ipt] + yp[ipt-1] ) / ( 0.5*( xp[ipt+1] - xp[ipt-1] ) )**2, \
    #                    ( yp[   -1] - 2*yp[ -2] + yp[   -3] ) / ( 0.5*( xp[   -1] - xp[   -3] ) )**2 ] )

    # ret = interpolate__polynomial( xp=xp, xd=xd, fx=fx, fpx=fpx, nDerivative=1, order=5 )

    
    # ------------------------------------------------- #
    # --- [3] case2: use ~2nd derivatives and 2 pts --- #
    # ------------------------------------------------- #
    xd   = np.array( [ xp[1], xp[-2] ] )
    fx   = np.array( [ yp[1], yp[-2] ] )
    fpx  = np.array( [ ( yp[    2] - yp[    0]) / ( xp[    2] - xp[    0] ), \
                       ( yp[   -1] - yp[   -3]) / ( xp[   -1] - xp[   -3] ) ] )
    fppx = np.array( [ ( yp[    2] - 2*yp[  1] + yp[    0] ) / ( 0.5*( xp[    2] - xp[    0] ) )**2, \
                       ( yp[   -1] - 2*yp[ -2] + yp[   -3] ) / ( 0.5*( xp[   -1] - xp[   -3] ) )**2 ] )

    ret = interpolate__polynomial( xp=xp, xd=xd, fx=fx, fpx=fpx, fppx=fppx, nDerivative=2, order=5 )

    # # ------------------------------------------------- #
    # # --- [4] case3: 3rd order= cubic interpolation --- #
    # # ------------------------------------------------- #
    # xd   = np.array( [ xp[1], xp[-2] ] )
    # fx   = np.array( [ yp[1], yp[-2] ] )
    # fpx  = np.array( [ ( yp[    2] - yp[    0]) / ( xp[    2] - xp[    0] ), \
    #                    ( yp[   -1] - yp[   -3]) / ( xp[   -1] - xp[   -3] ) ] )
    # fppx = np.array( [ ( yp[    2] - 2*yp[  1] + yp[    0] ) / ( 0.5*( xp[    2] - xp[    0] ) )**2, \
    #                    ( yp[   -1] - 2*yp[ -2] + yp[   -3] ) / ( 0.5*( xp[   -1] - xp[   -3] ) )**2 ] )

    # ret = interpolate__polynomial( xp=xp, xd=xd, fx=fx, fpx=fpx, nDerivative=1, order=3 )

    
    # ------------------------------------------------- #
    # --- [4] output result                         --- #
    # ------------------------------------------------- #
    diff = ret - yp
    stack = np.concatenate( [ret[:,None],yp[:,None],diff[:,None],], axis=1 )

    import nkUtilities.save__pointFile as spf
    outFile   = "test/out.dat"
    spf.save__pointFile( outFile=outFile, Data=stack )

    
    # ------------------------------------------------- #
    # --- [5] save in a png file                    --- #
    # ------------------------------------------------- #
    import nkUtilities.plot1D       as pl1
    import nkUtilities.load__config as lcf
    x_,y_                    = 0, 1
    pngFile                  = "test/out.png"
    config                   = lcf.load__config()
    config["plt_xAutoRange"] = True
    config["plt_yAutoRange"] = True
    config["plt_xRange"]     = [ -1.2, +1.2 ]
    config["plt_yRange"]     = [ -1.2, +1.2 ]

    fig     = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=xp, yAxis=ret, label="interpolated", linestyle="-" , linewidth=0.8 )
    fig.add__plot( xAxis=xp, yAxis=yp,  label="original"    , linestyle="--", linewidth=0.8 )
    fig.add__legend()
    fig.set__axis()
    fig.save__figure()


    

# ------------------------------------------------- #
# --- obsolete primitive ver.                   --- #
# ------------------------------------------------- #
# #  -- [1-1] 0th order                           --  #
# iDeriv = 0
# i_flr  = nCoef//(nDerivative+1) * iDeriv
# for ik in range( 0, nCoef//(nDerivative+1) ):
#     for jk in range( 0, nCoef ):
#         Amat[i_flr+ik,jk] = xd[ik]**( float(jk) )

# #  -- [1-2] 1st derivative order                --  #
# iDeriv = 1
# i_flr  = nCoef//(nDerivative+1) * iDeriv
# for ik in range( 0, nCoef//(nDerivative+1) ):
#     for jk in range( 1, nCoef ):
#         Amat[i_flr+ik,jk] = float(jk) * xd[ik]**( float(jk-1) )

# #  -- [1-3] 2nd derivative order                --  #
# iDeriv = 2
# i_flr  = nCoef//(nDerivative+1) * iDeriv
# for ik in range( 0, nCoef//(nDerivative+1) ):
#     for jk in range( 2, nCoef ):
#         Amat[i_flr+ik,jk] = float(jk) * float(jk-1) * xd[ik]**( float(jk-2) )


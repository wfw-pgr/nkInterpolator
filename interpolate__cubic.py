import numpy as np

# ----------------------------------------------- #
# --- cubic interpolation from end point      --- #
# ---  give (y1,y2) & (y'1, y'2) @ (x1,x2)    --- #
# ----------------------------------------------- #

# ========================================================= #
# ===  interpolate__cubic                               === #
# ========================================================= #

def interpolate__cubic( xp, x1, x2, f1, f2, fp1, fp2, extrapolate="cubic" ):

    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    if ( xp  is None ): sys.exit( "[interpolate__cubic] xp  == ???" )
    if ( x1  is None ): sys.exit( "[interpolate__cubic] x1  == ???" )
    if ( x2  is None ): sys.exit( "[interpolate__cubic] x2  == ???" )
    if ( f1  is None ): sys.exit( "[interpolate__cubic] f1  == ???" )
    if ( f2  is None ): sys.exit( "[interpolate__cubic] f2  == ???" )
    if ( fp1 is None ): sys.exit( "[interpolate__cubic] fp1 == ???" )
    if ( fp2 is None ): sys.exit( "[interpolate__cubic] fp2 == ???" )

    # ------------------------------------------------- #
    # --- [2] vector / Matrix Making                --- #
    # ------------------------------------------------- #
    #  -- [2-1] vectors                             --  #
    xv   = np.array( [ x1, x2 ] )
    lhs  = np.array( [ f1,f2,fp1,fp2] )
    #  -- [2-2] Amat                                --  #
    Amat = np.zeros( (4,4) )
    for i in range( 0, 2 ):
        for j in range( 0, 4 ):
            Amat[i,j] = xv[i]**( float(j) )
    for i in range( 2, 4 ):
        for j in range( 1, 4 ):
            Amat[i,j] = float(j) * xv[i-2]**( float(j-1) )

    # ------------------------------------------------- #
    # --- [3] get coefficients                      --- #
    # ------------------------------------------------- #
    coef = np.dot( np.linalg.inv( Amat ), lhs )

    # ------------------------------------------------- #
    # --- [4] interpolation                         --- #
    # ------------------------------------------------- #
    ret  = coef[3]*xp[:]**3 + coef[2]*xp[:]**2 + coef[1]*xp[:] + coef[0]

    
    # ------------------------------------------------- #
    # --- [5] linear extrapolation mode             --- #
    # ------------------------------------------------- #
    if ( extrapolate.lower() == "linear" ):
        ret = extrapolate__linear( xp, ret, x1, x2, f1, f2, fp1, fp2,  )

    return( ret )


# ========================================================= #
# ===  extrapolate__linear                              === #
# ========================================================= #
def extrapolate__linear( xp, yp, x1, x2, f1, f2, fp1, fp2 ):

    # ------------------------------------------------- #
    # --- [1] region determination                  --- #
    # ------------------------------------------------- #
    index1 = np.where( xp < x1 )
    index2 = np.where( xp > x2 )
    
    # ------------------------------------------------- #
    # --- [2] interpolation & extrapolation         --- #
    # ------------------------------------------------- #
    if ( len( index1[0] ) > 0 ):
        yp[index1] = fp1 * ( xp[index1] - x1 ) + f1
    if ( len( index2[0] ) > 0 ):
        yp[index2] = fp2 * ( xp[index2] - x2 ) + f2
    return( yp )

    

# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):

    # ------------------------------------------------- #
    # --- [1] Function parameter  ( End point )     --- #
    # ------------------------------------------------- #
    x1  = 0.0
    x2  = 1.0
    f1  = 0.0
    f2  = 1.0
    fp1 = 0.0
    fp2 = 2.0

    # ------------------------------------------------- #
    # --- [2] interpolation and analytical          --- #
    # ------------------------------------------------- #
    xp  = np.linspace( -1.0, 2.0, 31 )
    fp  = interpolate__cubic( xp=xp, x1=x1, x2=x2, f1=f1, f2=f2, fp1=fp1, fp2=fp2, \
                              extrapolate="linear" )
    xa  = np.linspace( -1.0, 2.0, 101 )
    fa  = xa**2

    # ------------------------------------------------- #
    # --- [3] plotting                              --- #
    # ------------------------------------------------- #
    import nkUtilities.plot1D       as pl1
    import nkUtilities.load__config as lcf
    pngFile = "png/out.png"
    config  = lcf.load__config()
    config["plt_xAutoRange"] = False
    config["plt_yAutoRange"] = False
    config["plt_xRange"]     = [-2.0, 3.0]
    config["plt_yRange"]     = [-1.0, 5.0]
    fig     = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=xa, yAxis=fa, label="analytical"   , linewidth=1.5 )
    fig.add__plot( xAxis=xp, yAxis=fp, label="interpolation", linewidth=0.0, \
                   marker="o", markersize=1.8 )
    fig.add__legend()
    fig.set__axis()
    fig.save__figure()

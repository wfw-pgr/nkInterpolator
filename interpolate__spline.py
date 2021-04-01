import numpy             as np
import scipy.interpolate as itp

# ---------------------------------------- #
# --- wrapper for spline interpolation --- #
# ---------------------------------------- #

# ========================================================= #
# ===  interpolate__spline.py                           === #
# ========================================================= #

def interpolate__spline( xval=None, xref=None, yref=None, kind="cubic" ):

    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    if ( xval is None ): sys.exit( "[interpolate__spline] xval == ???" )
    if ( xref is None ): sys.exit( "[interpolate__spline] xref == ???" )
    if ( yref is None ): sys.exit( "[interpolate__spline] yref == ???" )

    # ------------------------------------------------- #
    # --- [2] region determination                  --- #
    # ------------------------------------------------- #
    index1 = np.where(   xval <  xref[ 0] )
    index2 = np.where( ( xval >= xref[ 0] ) & ( xval <= xref[-1] ) )
    index3 = np.where(   xval >  xref[-1] )
    yval   = np.zeros_like( xval )
    
    # ------------------------------------------------- #
    # --- [3] interpolation & extrapolation         --- #
    # ------------------------------------------------- #
    if ( len( index1[0] ) > 0 ):
        dydx         = ( yref[1] - yref[0] ) / ( xref[1] - xref[0] )
        yval[index1] = dydx * ( xval[index1] - xref[0] ) + yref[0]
    if ( len( index3[0] ) > 0 ):
        dydx         = ( yref[-1] - yref[-2] ) / ( xref[-1] - xref[-2] )
        yval[index3] = dydx * ( xval[index3] - xref[-1] ) + yref[-1]
    if ( len( index2[0] ) > 0 ):
        interp_func  = itp.interp1d( xref, yref, kind=kind )
        yval[index2] = interp_func( xval[index2] )
    return( yval )


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):

    xref = np.linspace( 0.0, 3.0, 11 )
    yref = xref**2 + 2.0*xref
    
    xval = np.linspace( -1.0, 4.0, 31 )
    yval = interpolate__spline( xval=xval, xref=xref, yref=yref, kind="cubic" )

    xans = np.linspace( 0.0, 3.0, 101 )
    yans = xans**2 + 2.0*xans
    
    import nkUtilities.plot1D       as pl1
    import nkUtilities.load__config as lcf

    pngFile = "png/out.png"
    config  = lcf.load__config()
    fig     = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=xans, yAxis=yans, label="ans", linewidth=1.2 )
    fig.add__plot( xAxis=xref, yAxis=yref, label="ref", marker="o", markersize=1.8, linewidth=0.0 )
    fig.add__plot( xAxis=xval, yAxis=yval, label="val", marker="+", markersize=1.8, linewidth=0.0 )
    fig.add__legend()
    fig.set__axis()
    fig.save__figure()


    

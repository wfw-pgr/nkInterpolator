import os, sys
import numpy as np

# ========================================================= #
# ===  interpolate grided Data on line (2D)             === #
# ========================================================= #
def interpolate__gridData_onto_line( Data=None, gridFile="dat/grid.dat", lineFile=None, \
                                     x1  =[ 0.0, 0.0 ], x2=[ 1.0, 1.0 ], nDiv=101, \
                                     x_=0, y_=1, v_=2 ):

    # -- for vector field. -- #
    # x_,y_,v_ = 0, 1, 5
    # -- 
    
    # ------------------------------------------------- #
    # --- [1] grid & line Loading                   --- #
    # ------------------------------------------------- #
    #  -- [1-1] grid Data Loading                   --  #
    if   ( Data is not None ):
        gridData = np.copy( Data )
    elif ( gridFile is not None ):
        import nkUtilities.load__pointFile as lpf
        gridData  = lpf.load__pointFile( inpFile=gridFile, returnType="structured" )
        gridData  = gridData[:,:,(x_,y_,v_)]
    else:
        print( "[interpolate__gridData_onto_line] no Data & no gridFile... [ERROR]" )
        sys.exit()
    gridData_ = ( np.reshape( gridData, (-1,3) ) )
    
    #  -- [1-2] line Data Loading                   --  #
    if ( lineFile is None ):
        lineData  = generate__line2d( x1=x1, x2=x2, nDiv=nDiv )
    else:
        lineData  = lpf.load__pointFile( inpFile=lineFile, returnType="point" )
    lineData_        = np.zeros( (lineData.shape[0],3) )
    lineData_[:,0:2] = np.copy( lineData[:,0:2] )
    
    # ------------------------------------------------- #
    # --- [2] linear Interpolation                  --- #
    # ------------------------------------------------- #
    import nkInterpolator.LinearInterp2D as li2
    ret = li2.LinearInterp2D( gridData    =gridData, pointData=lineData_, \
                              gridDataType="structured" )
    
    # ------------------------------------------------- #
    # --- [3] grid Data contouring                  --- #
    # ------------------------------------------------- #
    import nkUtilities.plot1D     as pl1
    import nkUtilities.LoadConfig as lcf
    config = lcf.LoadConfig()
    config["plt_xAutoRange"] = True
    config["plt_yAutoRange"] = True
    config["plt_xRange"]     = [ 0.0,+1.0]
    config["plt_yRange"]     = [-1.0,+0.0]
    pl1.plot1D( xAxis=ret[:,0], yAxis=ret[:,2], pngFile="interpolate__gridData_onto_line.png", config=config )
    return( ret )
    

# ========================================================= #
# ===  sample line generator                            === #
# ========================================================= #

def generate__line2d( x1=[0.0,0.0], x2=[1.0,1.0], nDiv=101 ):

    ret       = np.zeros( (nDiv,2) )
    tval      = np.linspace( 0.0, 1.0, nDiv )
    ret[:,0]  = ( x2[0] - x1[0] ) * tval + x1[0]
    ret[:,1]  = ( x2[1] - x1[1] ) * tval + x1[1]
    return( ret )


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    import nkUtilities.generate__testprofile as gtp
    x1MinMaxNum = [ 0.0, 1.0, 11 ]
    x2MinMaxNum = [ 0.0, 1.0, 11 ]
    ret         = gtp.generate__testprofile( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
    	                                     returnType = "structured" )

    x1 = [0.0,0.0]
    x2 = [1.0,0.0]
    interpolate__gridData_onto_line( Data=ret, x1=x1, x2=x2 )

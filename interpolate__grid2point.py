import os, sys
import numpy as np


# ========================================================= #
# ===  interpolate__grid2point.py                       === #
# ========================================================= #

def interpolate__grid2point( gridData=None, pointData=None, dim=None, method="linear", size=None ):
    # ------------------------------------------------- #
    # --- [1] Argument                              --- #
    # ------------------------------------------------- #
    if ( gridData  is None ): sys.exit( "[interpolate__grid2point] gridData   == ???" )
    if ( pointData is None ): sys.exit( "[interpolate__grid2point] pointData  == ???" )

    if ( size is not None ):
        gridData = np.reshape( gridData, size )
    
    if ( dim is None ):
        if   ( gridData.ndim == 4 ):
            dim = 3
        elif ( gridData.ndim == 3 ):
            dim = 2
        else:
            print( "[interpolate__grid2point] unknown dimension :: dim = ???" )
            sys.exit()
    
    # ------------------------------------------------- #
    # --- [2]  call interpolation                   --- #
    # ------------------------------------------------- #

    if   ( dim == 2 ):

        if   ( method == "linear" ):
            import nkInterpolator.interpolate__bilinear as bil
            ret = bil.interpolate__bilinear( gridData=gridData, pointData=pointData )

        elif ( method == "cubic"  ):
            import nkInterpolator.interpolate__bicubic  as bic
            ret = bic.interpolate__bicubic ( gridData=gridData, pointData=pointData )

    elif ( dim == 3 ):

        if   ( method == "linear" ):
            import nkInterpolator.interpolate__trilinear as trl
            ret = trl.interpolate__trilinear( gridData=gridData, pointData=pointData )
            
        elif ( method == "cubic"  ):
            import nkInterpolator.interpolate__tricubic  as trc
            ret = trc.interpolate__tricubic ( gridData=gridData, pointData=pointData )


    return( ret )


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):

    method = "cubic"
    
    # ------------------------------------------------- #
    # --- [1] test ( 2D )                           --- #
    # ------------------------------------------------- #
    
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ 0.0, 1.0, 11 ]
    x2MinMaxNum = [ 0.0, 1.0, 11 ]
    x3MinMaxNum = [ 0.0, 0.0,  1 ]
    xRef        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    xRef[...,2] = np.sqrt( xRef[...,0]**2 + xRef[...,1]**2 )
    xRef        = np.reshape( xRef, (x2MinMaxNum[2],x1MinMaxNum[2],3) )
    
    x1MinMaxNum = [ 0.3, 0.7, 21 ]
    x2MinMaxNum = [ 0.3, 0.7, 21 ]
    x3MinMaxNum = [ 0.0, 0.0,  1 ]
    xItp        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )

    
    ret         = interpolate__grid2point( gridData=xRef, pointData=xItp, \
                                           method=method, dim=2 )

    Data        = np.zeros( (ret.shape[0],5) )
    Data[:,0:3] = ret[:,0:3]
    Data[:,  3] = np.sqrt( Data[:,0]**2 + Data[:,1]**2 )
    Data[:,  4] = Data[:,2] - Data[:,3]

    outFile     = "out_2d.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )

    print()
    print( " max( abs( diff ) ) == {0}".format( np.max( np.abs( Data[:,4] ) ) ) )
    print()

    

    # ------------------------------------------------- #
    # --- [2] test ( 3D )                           --- #
    # ------------------------------------------------- #
    
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ 0.0, 1.0, 11 ]
    x2MinMaxNum = [ 0.0, 1.0, 11 ]
    x3MinMaxNum = [ 0.0, 1.0, 11 ]
    xRef_       = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    xRef        = np.zeros( (x3MinMaxNum[2],x2MinMaxNum[2],x1MinMaxNum[2],4) )
    xRef[...,:3]= xRef_
    xRef[..., 3]= np.sqrt( xRef_[...,0]**2 + xRef_[...,1]**2 + xRef_[...,2]**2 )
    x1MinMaxNum = [ 0.3, 0.7, 21 ]
    x2MinMaxNum = [ 0.3, 0.7, 21 ]
    x3MinMaxNum = [ 0.3, 0.7, 21 ]
    xItp_       = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    xItp        = np.zeros( (x3MinMaxNum[2]*x2MinMaxNum[2]*x1MinMaxNum[2],4) )
    xItp[:,:3]  = np.copy( xItp_ )

    ret         = interpolate__grid2point( gridData=xRef, pointData=xItp, \
                                           method=method, dim=3 )
    
    Data        = np.zeros( (ret.shape[0],6) )
    Data[:,0:4] = ret[:,0:4]
    Data[:,  4] = np.sqrt( Data[:,0]**2 + Data[:,1]**2 + Data[:,2]**2 )
    Data[:,  5] = Data[:,3] - Data[:,4]

    outFile     = "out_3d.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )

    print()
    print( " max( abs( diff ) ) == {0}".format( np.max( np.abs( Data[:,5] ) ) ) )
    print()

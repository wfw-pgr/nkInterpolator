import numpy            as np
import numpy.ctypeslib  as Flib
import ctypes, sys
import os.path

# ================================================================ #
# ===  Linear Interpolation (2D)                               === #
# ================================================================ #
def interpolate__linear2D( gridData    =None , pointData=None, size=None, gridDataType="structured", \
                           gridDataOnly=False, unstructuredOrder="ijk" ):
    # ---------------------------------------- #
    # --- [1]   引数チェック               --- #
    # ---------------------------------------- #
    if ( gridData  is None ): sys.exit( "[interpolate__linear2D] gridData  ???" )
    if ( pointData is None ): sys.exit( "[interpolate__linear2D] pointData ???" )
    
    # ------------------------------------------------- #
    # --- [2] dx, dy, LI, LJ check                  --- #
    # ------------------------------------------------- #
    x_, y_, v_ = 0, 1, 2
    if   ( gridDataType == "structured"   ):
        if ( gridDataOnly ):
            # gridData :: [LJ,LI]
            if ( ( dx is None ) or ( dy is None ) or ( size is None ) ):
                # -- if  grid info. is unknown.... [ERROR] -- #
                print( "[interpolate__linear2D] gridDataType == unstructured, with gridDataOnly == True" )
                print( "[interpolate__linear2D]           dx == {0}".format(   dx ) )
                print( "[interpolate__linear2D]           dy == {0}".format(   dy ) )
                print( "[interpolate__linear2D]         size == {0}".format( size ) )
                sys.exit()
            else:
                # -- if  grid info. is given [OK]          -- #
                LI, LJ    = size[1], size[0]
                gridData_ = np.copy( gridData )
        else:
            # gridData :: [LJ,LI,3] :: (LJ,LI) * ( x_, y_, v_ )
            LJ, LI     = gridData.shape[0], gridData.shape[1]
            xMin, yMin = gridData[0,0,x_], gridData[0,0,y_]
            dx  , dy   = gridData[0,1,x_] - gridData[0,0,x_], gridData[1,0,y_] - gridData[0,0,y_]
            gridData_  = np.copy( gridData[:,:,v_] )

    elif ( gridDataType == "unstructured" ):
        # gridData :: [nData,3] = [LI*LJ,3]
        if ( size is None ):
            print( "[interpolate__linear2D] gridDataType == unstructured, but no size !!" )
            print( "[interpolate__linear2D]         size == ???" )
            sys.exit()
        else:
            LI, LJ   = size[0], size[1]
            if   ( unstructuredOrder == "ijk" ):
                # assume ascending ijk-order
                xMin  = gridData[ 0,x_]
                yMin  = gridData[ 0,y_]
                dx    = gridData[ 1,x_] - gridData[ 0,x_]
                dy    = gridData[LI,y_] - gridData[ 0,y_]
            elif ( unstructuredOrder == "kji" ):
                # assume ascending kji-order
                xMin  = gridData[ 0,x_]
                yMin  = gridData[ 0,y_]
                dx    = gridData[LJ,x_] - gridData[ 0,x_]
                dy    = gridData[ 1,y_] - gridData[ 0,y_]
            gridData_ = np.copy( gridData[:,v_].reshape( (LJ,LI) ) )
    nData = pointData.shape[0]
    
    # ---------------------------------------- #
    # --- [2]   引数準備                   --- #
    # ---------------------------------------- #
    #  -- [2-1] 使用する引数を準備         --  #
    #  -- [2-2] Fortranサイズへ変換        --  #
    gridData_  = np.array( gridData_, dtype=np.float64  )
    pointData_ = np.array( pointData, dtype=np.float64  )
    xMin_      = np.array(      xMin, dtype=np.float64  )
    yMin_      = np.array(      yMin, dtype=np.float64  )
    dx_        = np.array(      dx  , dtype=np.float64  )
    dy_        = np.array(      dy  , dtype=np.float64  )
    LI_        = ctypes.byref( ctypes.c_int64( LI )     )
    LJ_        = ctypes.byref( ctypes.c_int64( LJ )     )
    nData_     = ctypes.byref( ctypes.c_int64( nData )  )

    # ---------------------------------------- #
    # --- [3]   ライブラリをロード         --- #
    # ---------------------------------------- #
    #  -- [3-1] ライブラリを定義           --  #
    path   = os.path.expanduser('~') + "/.python/lib/nkInterpolator"
    pyLIB  = Flib.load_library( 'pylib.so', path )
    #  -- [3-2] 入出力管理                 --  #
    pyLIB.linearinterp2d_.argtypes = [
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.linearinterp2d_.restype = ctypes.c_void_p

    # ---------------------------------------- #
    # --- [4]   関数呼出 / 返却            --- #
    # ---------------------------------------- #
    pyLIB.linearinterp2d_( gridData_, pointData_, xMin_, yMin_, dx_, dy_, LI_, LJ_, nData_ )
    return( pointData_ )


# ================================================================ #
# ===  テスト用 呼び出し                                       === #
# ================================================================ #
if ( __name__=='__main__' ):
    import nkUtilities.equiSpaceGrid2D as esg
    
    # -- generate gridData  -- #
    x1MinMaxNum = [ -1.0, 1.0, 11 ]
    x2MinMaxNum = [ -1.0, 1.0, 11 ]
    grid = esg.equiSpaceGrid2D( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                returnType="point" )
    gridData      = np.zeros( (grid.shape[0],3) )
    gridData[:,0] = grid[:,0]
    gridData[:,1] = grid[:,1]
    gridData[:,2] = 1.0 - ( grid[:,0]**2 + grid[:,1]**2 )

    import nkUtilities.cMapTri as cmt
    cmt.cMapTri( xAxis=gridData[:,0], yAxis=gridData[:,1], cMap=gridData[:,2], \
                 pngFile="xRef.png" )

    # -- generate pointData -- #
    x1MinMaxNum = [ -1.0, 1.0, 101 ]
    x2MinMaxNum = [ -1.0, 1.0, 51 ]
    point = esg.equiSpaceGrid2D( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                 returnType="point" )
    pointData      = np.zeros( (point.shape[0],3) )
    pointData[:,0] = point[:,0]
    pointData[:,1] = point[:,1]
    np.savetxt( "tmp.dat", pointData )
    print( pointData.shape )

    ret = interpolate__linear2D( gridData=gridData, pointData=pointData, \
                                 size=(11,11), gridDataType="unstructured" )
    xAxis = ret[:,0]
    yAxis = ret[:,1]
    zAxis = ret[:,2]
    import nkUtilities.cMapTri as cmt
    cmt.cMapTri( xAxis=xAxis, yAxis=yAxis, cMap=zAxis, \
                 pngFile="xItp.png" )
    print( ret.shape )

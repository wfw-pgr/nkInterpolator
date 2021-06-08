import numpy            as np
import numpy.ctypeslib  as Flib
import ctypes, sys
import os.path

# ================================================================ #
# ===  interpolate__trilinear                                  === #
# ================================================================ #
def interpolate__trilinear( gridData=None, pointData=None, outOfRangeMode="zero" ):

    x_, y_, z_, v_ = 0, 1, 2, 3
    
    # ------------------------------------------------- #
    # --- [1]   引数チェック                        --- #
    # ------------------------------------------------- #
    if ( gridData  is None ): sys.exit( "[interpolate__trilinear] gridData  ???" )
    if ( pointData is None ): sys.exit( "[interpolate__trilinear] pointData ???" )

    # ------------------------------------------------- #
    # --- [2]   引数準備                            --- #
    # ------------------------------------------------- #
    #  -- [2-1] 使用する引数を準備                  --  #
    gData      = np.copy( gridData[:,:,:,v_] )
    pData      = np.copy( pointData )
    LI,LJ,LK   = gridData.shape[2],gridData.shape[1],gridData.shape[0]
    nData      = pointData.shape[0]
    xMin       = gridData[0,0,0,x_]
    yMin       = gridData[0,0,0,y_]
    zMin       = gridData[0,0,0,z_]
    dx         = ( gridData[0,0,1,x_] - gridData[0,0,0,x_] )
    dy         = ( gridData[0,1,0,y_] - gridData[0,0,0,y_] )
    dz         = ( gridData[1,0,0,z_] - gridData[0,0,0,z_] )
    if   ( outOfRangeMode == "forced" ):
        i_outOfRangeMode = 1
    elif ( outOfRangeMode == "zero"   ):
        i_outOfRangeMode = 0
    elif ( outOfRangeMode == "large"  ):
        i_outOfRangeMode = 2
    else:
        print( "[interpolate__trilinear.py] outOfRangeMode == {0} ??? [ foreced / zero / large ]".format( outOfRangeMode ) )
        sys.exit()
    
    #  -- [2-2] Fortranサイズへ変換                 --  #
    gData_     =     np.array( gData      , dtype=np.float64  )
    pData_     =     np.array( pData      , dtype=np.float64  )
    xMin_      =     np.array( xMin       , dtype=np.float64  )
    yMin_      =     np.array( yMin       , dtype=np.float64  )
    zMin_      =     np.array( zMin       , dtype=np.float64  )
    dx_        =     np.array( dx         , dtype=np.float64  )
    dy_        =     np.array( dy         , dtype=np.float64  )
    dz_        =     np.array( dz         , dtype=np.float64  )
    LI_        = ctypes.byref( ctypes.c_int64( LI   )  )
    LJ_        = ctypes.byref( ctypes.c_int64( LJ   )  )
    LK_        = ctypes.byref( ctypes.c_int64( LK   )  )
    nData_     = ctypes.byref( ctypes.c_int64( nData )  )
    i_outOfRangeMode_  = ctypes.byref( ctypes.c_int64( i_outOfRangeMode )  )
    
    # ------------------------------------------------- #
    # --- [3]   ライブラリをロード                  --- #
    # ------------------------------------------------- #
    
    #  -- [3-1] ライブラリを定義                    --  #
    pyLIB  = Flib.load_library( 'pylib.so', os.path.abspath( os.path.dirname(__file__) ) )
    
    #  -- [3-2] 入出力管理                          --  #
    pyLIB.interpolate__trilinear_.argtypes = [
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.interpolate__trilinear_.restype = ctypes.c_void_p

    # ------------------------------------------------- #
    # --- [4]   関数呼出 / 返却                     --- #
    # ------------------------------------------------- #
    pyLIB.interpolate__trilinear_( gData_, pData_, xMin_, yMin_, zMin_, \
                                   dx_, dy_, dz_, LI_, LJ_, LK_, nData_, i_outOfRangeMode_ )
    return( pData_ )


# ================================================================ #
# ===  テスト用 呼び出し                                       === #
# ================================================================ #
if ( __name__=='__main__' ):

    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ 0.0, 1.0, 11 ]
    x2MinMaxNum = [ 0.0, 1.0, 11 ]
    x3MinMaxNum = [ 0.0, 1.0, 11 ]
    gData_      = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                      x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    gData       = np.zeros( (x3MinMaxNum[2],x2MinMaxNum[2],x1MinMaxNum[2],4) )
    gData[...,:3]= gData_
    gData[..., 3]= np.sqrt( gData_[...,0]**2 + gData_[...,1]**2 + gData_[...,2]**2 )
    x1MinMaxNum = [ 0.3, 0.7, 21 ]
    x2MinMaxNum = [ 0.3, 0.7, 21 ]
    x3MinMaxNum = [ 0.3, 0.7, 21 ]
    pData_      = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                      x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    pData       = np.zeros( (x3MinMaxNum[2]*x2MinMaxNum[2]*x1MinMaxNum[2],4) )
    pData[:,:3] = np.copy( pData_ )

    ret         = interpolate__trilinear( gridData=gData, pointData=pData )
    
    Data        = np.zeros( (ret.shape[0],5) )
    Data[:,0:4] = ret[:,0:4]
    Data[:,  4] = np.sqrt( Data[:,0]**2 + Data[:,1]**2 + Data[:,2]**2 )

    outFile     = "out.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )


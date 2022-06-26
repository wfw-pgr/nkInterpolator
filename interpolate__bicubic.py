import numpy            as np
import numpy.ctypeslib  as Flib
import ctypes, sys
import os.path

# ------------------------------------------------ #
# --   * gridData  :: [LJ,LI,3]      (x,y,v)    -- #
# --   * pointData :: [npt,3]                   -- #
# ------------------------------------------------ #


# ================================================================ #
# ===  interpolate__bicubic                                    === #
# ================================================================ #
def interpolate__bicubic( gridData=None, pointData=None, LK__no_of_copy_layer=5 ):

    x_,y_,z_,v_ = 0, 1, 2, 3
    
    # ------------------------------------------------- #
    # --- [1]   引数チェック                        --- #
    # ------------------------------------------------- #
    if ( gridData  is None ): sys.exit( "[interpolate__bicubic] gridData  ???" )
    if ( pointData is None ): sys.exit( "[interpolate__bicubic] pointData ???" )

    # ------------------------------------------------- #
    # --- [2]   引数準備                            --- #
    # ------------------------------------------------- #
    #  -- [2-1] 使用する引数を準備                  --  #
    LI,LJ,LK   = gridData.shape[1],gridData.shape[0], LK__no_of_copy_layer
    nItp       = pointData.shape[0]

    xRef_      = np.zeros( (LK,LJ,LI,4) )
    xItp_      = np.zeros( (nItp,4)     )
    zcoord     = np.linspace( -2.0, 2.0, LK__no_of_copy_layer )

    for ik in range( LK ):
        xRef_[ik,:,:,x_] = np.copy( gridData[:,:,0] )
        xRef_[ik,:,:,y_] = np.copy( gridData[:,:,1] )
        xRef_[ik,:,:,z_] = zcoord[ik]
        xRef_[ik,:,:,v_] = np.copy( gridData[:,:,2] )
        xItp_[:,x_]      = np.copy( pointData[:,0] )
        xItp_[:,y_]      = np.copy( pointData[:,1] )
        xItp_[:,z_]      = 0.0

    #  -- [2-2] Fortranサイズへ変換                 --  #
    xRef_      =     np.array( xRef_ , dtype=np.float64  )
    xItp_      =     np.array( xItp_ , dtype=np.float64  )
    LI_        = ctypes.byref( ctypes.c_int64( LI   )  )
    LJ_        = ctypes.byref( ctypes.c_int64( LJ   )  )
    LK_        = ctypes.byref( ctypes.c_int64( LK   )  )
    nItp_      = ctypes.byref( ctypes.c_int64( nItp )  )
    
    # ------------------------------------------------- #
    # --- [3]   ライブラリをロード                  --- #
    # ------------------------------------------------- #
    
    #  -- [3-1] ライブラリを定義                    --  #
    pyLIB  = Flib.load_library( 'pylib.so', os.path.abspath( os.path.dirname(__file__) ) )
    
    #  -- [3-2] 入出力管理                          --  #
    pyLIB.cubicinterpolation_3d_.argtypes = [
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.cubicinterpolation_3d_.restype = ctypes.c_void_p

    # ------------------------------------------------- #
    # --- [4]   関数呼出 / 返却                     --- #
    # ------------------------------------------------- #
    pyLIB.cubicinterpolation_3d_( xRef_, xItp_, LI_,LJ_,LK_,nItp_, )
    ret      = np.zeros( (nItp,3) )
    ret[:,0] = np.copy ( xItp_[:,0] )
    ret[:,1] = np.copy ( xItp_[:,1] )
    ret[:,2] = np.copy ( xItp_[:,3] )
    return( ret )


# ================================================================ #
# ===  テスト用 呼び出し                                       === #
# ================================================================ #
if ( __name__=='__main__' ):

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
    ret         = interpolate__bicubic( gridData=xRef, pointData=xItp )
    Data        = np.zeros( (ret.shape[0],5) )
    Data[:,0:3] = np.copy( ret[:,0:3] )
    Data[:,  3] = np.sqrt( Data[:,0]**2 + Data[:,1]**2 )
    Data[:,  4] = Data[:,2] - Data[:,3]

    outFile     = "out.dat"
    names       = [ "x", "y", "interped", "answer", "diff" ]
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )


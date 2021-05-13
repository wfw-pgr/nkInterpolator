import numpy            as np
import numpy.ctypeslib  as Flib
import ctypes, sys
import os.path

# ================================================================ #
# ===  interpolate__tricubic                                   === #
# ================================================================ #
def interpolate__tricubic( xRef=None, xItp=None ):
    # ------------------------------------------------- #
    # --- [1]   引数チェック                        --- #
    # ------------------------------------------------- #
    if ( xRef is None ): sys.exit( "[interpolate__tricubic] xRef ???" )
    if ( xItp is None ): sys.exit( "[interpolate__tricubic] xItp ???" )

    # ------------------------------------------------- #
    # --- [2]   引数準備                            --- #
    # ------------------------------------------------- #
    #  -- [2-1] 使用する引数を準備                  --  #
    LI,LJ,LK = xRef.shape[2],xRef.shape[1],xRef.shape[0]
    nItp     = xItp.shape[0]
    
    #  -- [2-2] Fortranサイズへ変換                 --  #
    xRef_      =     np.array( xRef  , dtype=np.float64  )
    xItp_      =     np.array( xItp  , dtype=np.float64  )
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
    return( xItp_ )


# ================================================================ #
# ===  テスト用 呼び出し                                       === #
# ================================================================ #
if ( __name__=='__main__' ):

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

    ret         = interpolate__tricubic( xRef=xRef, xItp=xItp )
    
    Data        = np.zeros( (ret.shape[0],5) )
    Data[:,0:4] = ret[:,0:4]
    Data[:,  4] = np.sqrt( Data[:,0]**2 + Data[:,1]**2 + Data[:,2]**2 )

    outFile     = "out.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )


import numpy            as np
import numpy.ctypeslib  as Flib
import ctypes, sys
import os.path


# ================================================================ #
# ===  linear Interpolation 1D ( non-uniform grid ver. )       === #
# ================================================================ #
def LinearInterp1D( xa=None, fa=None, xp=None, force_to_interpolate=False ):
    # ---------------------------------------- #
    # --- [1]   引数チェック               --- #
    # ---------------------------------------- #
    if ( xa is None ): sys.exit( "[linearInterp1D_nonUniform] xa ???" )
    if ( fa is None ): sys.exit( "[linearInterp1D_nonUniform] fa ???" )
    if ( xp is None ): sys.exit( "[linearInterp1D_nonUniform] xp ???" )
    
    # ---------------------------------------- #
    # --- [2]   引数準備                   --- #
    # ---------------------------------------- #
    #  -- [2-1] 使用する引数を準備         --  #
    nData    = xa.shape[0]
    nIntp    = xp.shape[0]
    fp       = np.zeros( ( nIntp, ) )
    Flag     = int( force_to_interpolate )
    #  -- [2-2] Fortranサイズへ変換        --  #
    xa_      =     np.array( xa  , dtype=np.float64  )
    fa_      =     np.array( fa  , dtype=np.float64  )
    xp_      =     np.array( xp  , dtype=np.float64  )
    fp_      =     np.array( fp  , dtype=np.float64  )
    nData_   = ctypes.byref( ctypes.c_int64( nData ) )
    nIntp_   = ctypes.byref( ctypes.c_int64( nIntp ) )
    Flag_    = ctypes.byref( ctypes.c_int64( Flag  ) )

    # ---------------------------------------- #
    # --- [3]   ライブラリをロード         --- #
    # ---------------------------------------- #
    #  -- [3-1] ライブラリを定義           --  #
    path   = os.path.expanduser('~') + "/.python/lib/interpRoutines"
    pyLIB  = Flib.load_library( 'pylib.so', path )
    #  -- [3-2] 入出力管理                 --  #
    pyLIB.linearinterp1d_.argtypes = [
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.linearinterp1d_.restype = ctypes.c_void_p

    # ---------------------------------------- #
    # --- [4]   関数呼出 / 返却            --- #
    # ---------------------------------------- #
    pyLIB.linearinterp1d_( xa_, fa_, xp_, fp_, nData_, nIntp_, Flag_ )
    return( fp_ )


# ================================================================ #
# ===  テスト用 呼び出し                                       === #
# ================================================================ #
if ( __name__=='__main__' ):

    xa = np.linspace( 0.0, 2.0,  11 )
    xp = np.linspace( 0.0, 2.0, 101 )
    fa = xa**2
    fp = LinearInterp1D( xa=xa, xp=xp, fa=fa )

    import nkUtilities.plot1D as pl1
    fig = pl1.plot1D( pngFile="png/out.png" )
    fig.add__plot( xAxis=xa, yAxis=fa, label="fa" )
    fig.add__plot( xAxis=xp, yAxis=fp, label="fp", linestyle="--" )
    fig.set__axis()
    fig.save__figure()

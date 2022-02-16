import numpy           as np
import numpy.ctypeslib as clib
import ctypes, sys
import os.path

# ================================================================ #
# ===  interpolate__fromTriMesh                                === #
# ================================================================ #
def interpolate__fromTriMesh( elems=None, nodes=None, points=None, index_begin_from_zero=True ):
    
    # ------------------------------------------------- #
    # --- [1]   引数チェック                        --- #
    # ------------------------------------------------- #
    if ( elems   is None ): sys.exit( "[interpolate__fromTriMesh] elems   ???" )
    if ( nodes   is None ): sys.exit( "[interpolate__fromTriMesh] nodes   ???" )
    if ( points  is None ): sys.exit( "[interpolate__fromTriMesh] points  ???" )

    if ( index_begin_from_zero ):
        elems[:,:] = elems[:,:] + 1
    
    # ------------------------------------------------- #
    # --- [2]   引数準備                            --- #
    # ------------------------------------------------- #
    #  -- [2-1] 使用する引数を準備                  --  #
    nElems     =  elems.shape[0]
    nNodes     =  nodes.shape[0]
    nPoints    = points.shape[0]
    #  -- [2-2] Fortranサイズへ変換                 --  #
    nodes_     =     np.array( nodes   , dtype=np.float64 )
    elems_     =     np.array( elems   , dtype=np.int32   )
    points_    =     np.array( points  , dtype=np.float64 )
    nNodes_    = ctypes.byref( ctypes.c_int64( nNodes  )  )
    nElems_    = ctypes.byref( ctypes.c_int64( nElems  )  )
    nPoints_   = ctypes.byref( ctypes.c_int64( nPoints )  )
    
    # ------------------------------------------------- #
    # --- [3]   ライブラリをロード                  --- #
    # ------------------------------------------------- #
    #  -- [3-1] ライブラリを定義                    --  #
    pyLIB  = clib.load_library( 'pylib.so', os.path.abspath( os.path.dirname(__file__) ) )
    
    #  -- [3-2] 入出力管理                          --  #
    pyLIB.interpolate__fromtrimesh_.argtypes = [
        clib.ndpointer( dtype=np.float64 ),
        clib.ndpointer( dtype=np.int32   ),
        clib.ndpointer( dtype=np.float64 ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.interpolate__fromtrimesh_.restype = ctypes.c_void_p

    # ------------------------------------------------- #
    # --- [4]   関数呼出 / 返却                     --- #
    # ------------------------------------------------- #
    pyLIB.interpolate__fromtrimesh_( nodes_, elems_, points_, nNodes_, nElems_, nPoints_ )
    return( points_ )


# ================================================================ #
# ===  テスト用 呼び出し                                       === #
# ================================================================ #
if ( __name__=='__main__' ):

    x_, y_, z_ = 0, 1, 2

    elemFile  = "test/elems.dat"
    nodeFile  = "test/nodes.dat"
    pointFile = "test/points.dat"

    import nkUtilities.load__pointFile as lpf
    elems     = lpf.load__pointFile( inpFile=elemFile , returnType="point" )
    nodes     = lpf.load__pointFile( inpFile=nodeFile , returnType="point" )
    points    = lpf.load__pointFile( inpFile=pointFile, returnType="point" )

    nodes[:,z_] = np.sqrt( nodes[:,x_]**2 + nodes[:,y_]**2 )
    
    ret       = interpolate__fromTriMesh( elems=elems, nodes=nodes, points=points )
    print( ret )

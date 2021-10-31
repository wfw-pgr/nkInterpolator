import numpy            as np
import numpy.ctypeslib  as Flib
import ctypes, sys
import os.path
import scipy.spatial    as spa

# ================================================================ #
# ===  barycentric__interpolator                               === #
# ================================================================ #
def barycentric__interpolator( nodes=None, points=None ):
    # ---------------------------------------- #
    # --- [1]   引数チェック               --- #
    # ---------------------------------------- #
    if ( nodes  is None ): sys.exit( "[barycentric__interpolator] nodes  ???" )
    if ( points is None ): sys.exit( "[barycentric__interpolator] points ???" )

    # ---------------------------------------- #
    # --- [2]   引数準備                   --- #
    # ---------------------------------------- #
    #  -- [2-1] 使用する引数を準備         --  #

    delaunay  = spa.Delaunay( nodes[:,0:2] )
    simplex   = np.array( delaunay.simplices                    , dtype=np.int32 )
    pwhere    = np.array( delaunay.find_simplex( points[:,0:2] ), dtype=np.int32 )
    nNodes    =   nodes.shape[0]
    nPoints   =  points.shape[0]
    nSimplex  = simplex.shape[0]

    #  -- [2-2] Fortranサイズへ変換        --  #
    nodes_    =     np.array( nodes  , dtype=np.float64  )
    points_   =     np.array( points , dtype=np.float64  )
    simplex_  =     np.array( simplex, dtype=np.int32    )
    pwhere_   =     np.array( pwhere , dtype=np.int32    )
    nNodes_   = ctypes.byref( ctypes.c_int64( nNodes   ) )
    nPoints_  = ctypes.byref( ctypes.c_int64( nPoints  ) )
    nSimplex_ = ctypes.byref( ctypes.c_int64( nSimplex ) )

    # ---------------------------------------- #
    # --- [3]   ライブラリをロード         --- #
    # ---------------------------------------- #
    #  -- [3-1] ライブラリを定義           --  #
    path   = os.path.expanduser('~') + "/.python/lib/nkInterpolator"
    pyLIB  = Flib.load_library( 'pylib.so', path )
    #  -- [3-2] 入出力管理                 --  #
    pyLIB.barycentric__interpolator_.argtypes = [
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.int32   ),
        Flib.ndpointer( dtype=np.int32   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.barycentric__interpolator_.restype = ctypes.c_void_p

    # ---------------------------------------- #
    # --- [4]   関数呼出 / 返却            --- #
    # ---------------------------------------- #
    pyLIB.barycentric__interpolator_( nodes_ , points_ , simplex_, pwhere_, \
                                      nNodes_, nPoints_, nSimplex_ )

    # ------------------------------------------------- #
    # --- [5] debug mode                            --- #
    # ------------------------------------------------- #
    # debugMode = False
    # if ( debugMode ):
    #     with open( "debug_pwhere.dat" , "w" ) as f:
    #         np.savetxt( f, pwhere )
    #     with open( "debug_simplex.dat", "w" ) as f:
    #         np.savetxt( f, simplex )
    
    return( points_ )



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):

    # ------------------------------------------------- #
    # --- [1] grid data making                      --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -1.0, 1.0, 201 ]
    x2MinMaxNum = [ -1.0, 1.0, 201 ]
    grid        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     returnType = "point" )
    radii       = np.sqrt( grid[:,0]**2 + grid[:,1]**2 )
    index       = np.where( radii <= 1.0 )
    grid        = grid [index]
    radii       = radii[index]
    height      = ( np.cos( 0.5*np.pi*radii ) )**2
    gData       = np.concatenate( [grid,np.reshape( height, (-1,1))], 1 )

    # ------------------------------------------------- #
    # --- [2] point data making                     --- #
    # ------------------------------------------------- #
    nPoints     = 10001
    points      = np.zeros( (nPoints,3) )
    points[:,0] = ( 1.0 - ( -1.0 ) ) * np.random.rand( nPoints ) + ( -1.0 ) 
    points[:,1] = ( 1.0 - ( -1.0 ) ) * np.random.rand( nPoints ) + ( -1.0 ) 
    radii       = np.sqrt( points[:,0]**2 + points[:,1]**2 )
    index       = np.where( radii < 1.0 )
    points      = points[index]

    # ------------------------------------------------- #
    # --- [3] interpolator                          --- #
    # ------------------------------------------------- #

    ret = barycentric__interpolator( nodes=gData, points=points )

    import nkUtilities.cMapTri      as cmt
    import nkUtilities.load__config as lcf
    config = lcf.load__config()
    config["cmp_AutoLevel"]  = False
    config["cmp_MaxMin"]     = [0.0,1.0]
    config["cmp_xAutoRange"] = False
    config["cmp_yAutoRange"] = False
    config["cmp_xRange"]     = [-1.0,+1.0]
    config["cmp_yRange"]     = [-1.0,+1.0]
    cmt.cMapTri( xAxis=gData[:,0], yAxis=gData[:,1], cMap=gData[:,2], pngFile="png/out1.png", config=config )
    cmt.cMapTri( xAxis=  ret[:,0], yAxis=  ret[:,1], cMap=  ret[:,2], pngFile="png/out2.png", config=config )

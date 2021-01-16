import numpy as np

# ========================================================= #
# ===  make boundary cell coordinates                   === #
# ========================================================= #


def make__boundary_cell():

    x_, y_, z_ = 0, 1, 2
    
    # ------------------------------------------------- #
    # --- [1] coordinate settings                   --- #
    # ------------------------------------------------- #

    import nkUtilities.load__constants as lcn
    cnsFile  = "dat/parameter.conf"
    const    = lcn.load__constants( inpFile=cnsFile )
    
    import nkUtilities.equiSpaceGrid as esg
    grid     = esg.equiSpaceGrid( x1MinMaxNum=const["x1MinMaxNum"], x2MinMaxNum=const["x2MinMaxNum"], \
                                  x3MinMaxNum=const["x3MinMaxNum"], returnType = "point" )
    BField        = np.zeros( (grid.shape[0],6) )
    BField[:,0:3] = np.copy( grid )

    # ------------------------------------------------- #
    # --- [2] boundary detection                    --- #
    # ------------------------------------------------- #

    nField           = BField.shape[0]
    LI,LJ,LK         = int( const["x1MinMaxNum"][2] ), int( const["x2MinMaxNum"][2] ), int( const["x3MinMaxNum"][2] )
    dx               = ( const["x1MinMaxNum"][1] - const["x1MinMaxNum"][0] ) / ( LI-1 )
    dy               = ( const["x2MinMaxNum"][1] - const["x2MinMaxNum"][0] ) / ( LJ-1 )
    delta            = np.min( [ dx, dy ] )
    radii            = np.sqrt( BField[:,x_]**2 + BField[:,y_]**2 )
    zpos             = np.copy( BField[:,z_] )
    r1, r2           = const["radius"], const["radius"] + delta
    z1, z2           = const["x3MinMaxNum"][0], const["x3MinMaxNum"][1]
    index1           = np.where( ( radii >= r1 ) & ( radii < r2 ) )
    index2           = np.where( ( radii <  r2 ) & ( zpos == z1 ) )
    index3           = np.where( ( radii <  r2 ) & ( zpos == z2 ) )
    boundary         = np.zeros( (nField,) )
    boundary[index1] = 1.0
    boundary[index2] = 1.0
    boundary[index3] = 1.0

    bdr_surface      = np.copy( boundary )

    index4           = np.where( ( radii >= r2 ) )
    boundary[index4] = 1.0
    bdr_flags        = np.copy( boundary )
    
    source           = np.zeros( (nField,8) )
    source[:,0:6]    = BField
    source[:,  6]    = bdr_surface
    source[:,  7]    = bdr_flags
    source           = np.reshape( source, (LK,LJ,LI,8) )

    index            = np.where( bdr_surface[:] == 1.0 )
    BField_boundary  = ( BField[index] ) [:,0:3]

    # ------------------------------------------------- #
    # --- [4] save in File                          --- #
    # ------------------------------------------------- #
    
    #  -- [4-1] save coordinate                     --  #
    outFile   = "dat/ems_pst.coord"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=BField_boundary )

    #  -- [4-2] save source                         --  #
    outFile   = "dat/source_wofield.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=source          )
    
    #  -- [4-3] display source                      --  #
    import nkVTKRoutines.convert__vtkStructuredGrid as vts
    outFile  = "png/source.vts"
    names    = [ "bx", "by", "bz", "boundary", "flag" ]
    vts.convert__vtkStructuredGrid( Data=source, outFile=outFile, names=names )



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    make__boundary_cell()

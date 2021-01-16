import os, sys
import numpy as np

# ========================================================= #
# ===  make__boundary.py                                === #
# ========================================================= #

def make__boundary():

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
    BField[:,0:3] = grid
    
    # ------------------------------------------------- #
    # --- [2] coil position settings                --- #
    # ------------------------------------------------- #

    x_,y_,z_ = 0, 1, 2
    theta    = np.linspace( 0.0, 2.0*np.pi, const["ntheta"] )
    coil_x   =   const["coil_center"][x_] + np.cos( theta )
    coil_y   =   const["coil_center"][y_] + np.sin( theta )
    coil_z1  = + const["coil_center"][z_] + theta * 0.0
    coil_z2  = - const["coil_center"][z_] + theta * 0.0
    coil1    = np.concatenate( [ coil_x[:,None], coil_y[:,None], coil_z1[:,None] ], axis=1 )
    coil2    = np.concatenate( [ coil_x[:,None], coil_y[:,None], coil_z2[:,None] ], axis=1 )
    
    # ------------------------------------------------- #
    # --- [3] coil field                            --- #
    # ------------------------------------------------- #

    import nkPhysicsRoutines.calc__biotSavartBField as bsf
    field1       = bsf.calc__biotSavartBField( bfield=BField, coils=coil1, I0=const["I0"] )
    field2       = bsf.calc__biotSavartBField( bfield=BField, coils=coil2, I0=const["I0"] )
    field        = np.zeros( (field1.shape[0],field1.shape[1]) )
    field[:,0:3] = field1[:,0:3]
    field[:,3:6] = field1[:,3:6] + field2[:,3:6]

    shape        = ( int( const["x3MinMaxNum"][2] ), int( const["x2MinMaxNum"][2] ), \
                     int( const["x1MinMaxNum"][2] ), 6 )
    field        = np.reshape( field, shape )
    
    # ------------------------------------------------- #
    # --- [4] save point Field                      --- #
    # ------------------------------------------------- #

    outFile      = "dat/coilfield.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=field )

    # ------------------------------------------------- #
    # --- [5] source with boundary info             --- #
    # ------------------------------------------------- #
    #  -- for laplacian ( R.H.S. == 0 for laplacian ) -- #
    LI, LJ, LK          = field.shape[2], field.shape[1], field.shape[0]
    source_             = np.copy( np.reshape( field, (-1,6) ) )
    radii               = np.sqrt( source_[:,x_]**2 + source_[:,y_]**2 )
    index               = np.where( radii <= const["radius"] )
    boundary            = np.ones( (source_.shape[0],) )
    boundary[index]     = 0.0
    source_[index,3:6]  = 0.0
    source_             = np.reshape( source_ , (LK,LJ,LI,6) )
    source_[ 0,:,:,3:6] = field[ 0,:,:,3:6]
    source_[-1,:,:,3:6] = field[-1,:,:,3:6]
    boundary            = np.reshape( boundary, (LK,LJ,LI,1) )
    
    srcFile      = "dat/source.dat"
    spf.save__pointFile( outFile=srcFile, Data=source_ )
    bdrFile      = "dat/boundary.dat"
    spf.save__pointFile( outFile=bdrFile, Data=boundary )
    
    # ------------------------------------------------- #
    # --- [7] convert vts File                      --- #
    # ------------------------------------------------- #
    import nkVTKRoutines.convert__vtkStructuredGrid as vts
    outFile  = "png/field.vts"
    names    = [ "bx", "by", "bz" ]
    vts.convert__vtkStructuredGrid( Data=field, outFile=outFile, names=names )

    outFile  = "png/source.vts"
    names    = [ "bx", "by", "bz" ]
    vts.convert__vtkStructuredGrid( Data=source_, outFile=outFile, names=names )

    return()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    make__boundary()

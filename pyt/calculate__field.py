import numpy as np


# ========================================================= #
# ===  calculate__field                                 === #
# ========================================================= #

def calculate__field():


    # ------------------------------------------------- #
    # --- [1] coordinate settings                   --- #
    # ------------------------------------------------- #

    import nkUtilities.load__constants as lcn
    cnsFile  = "dat/parameter.conf"
    const    = lcn.load__constants( inpFile=cnsFile )

    
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
    # --- [3] load BField coordinate                --- #
    # ------------------------------------------------- #

    inpFile       = "dat/ems_pst.coord"
    import nkUtilities.load__pointFile as lpf
    coord         = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    BField        = np.zeros( (coord.shape[0],6) )
    BField[:,0:3] = coord
    
    
    # ------------------------------------------------- #
    # --- [4] coil field                            --- #
    # ------------------------------------------------- #

    import nkPhysicsRoutines.calc__biotSavartBField as bsf
    field1       = bsf.calc__biotSavartBField( bfield=BField, coils=coil1, I0=const["I0"] )
    field2       = bsf.calc__biotSavartBField( bfield=BField, coils=coil2, I0=const["I0"] )
    field        = np.zeros( (field1.shape[0],field1.shape[1]) )
    field[:,0:3] = field1[:,0:3]
    field[:,3:6] = field1[:,3:6] + field2[:,3:6]

    # shape        = ( int( const["x3MinMaxNum"][2] ), int( const["x2MinMaxNum"][2] ), \
    #                  int( const["x1MinMaxNum"][2] ), 6 )
    # field        = np.reshape( field, shape )
    
    
    # ------------------------------------------------- #
    # --- [5] save point Field                      --- #
    # ------------------------------------------------- #

    outFile      = "dat/ems_pst.field"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=field )



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    calculate__field()

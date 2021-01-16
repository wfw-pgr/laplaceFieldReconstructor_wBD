import numpy as np

# ========================================================= #
# ===  Load BField from ems_pst.field onto source.dat   === #
# ========================================================= #

def load__field():

    # ------------------------------------------------- #
    # --- [1] load ems_pst.field & source.dat       --- #
    # ------------------------------------------------- #

    import nkUtilities.load__constants as lcn
    cnsFile  = "dat/parameter.conf"
    const    = lcn.load__constants( inpFile=cnsFile )

    import nkUtilities.load__pointFile as lpf
    inpFile = "dat/ems_pst.field"
    BField  = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    inpFile = "dat/source_wofield.dat"
    source  = lpf.load__pointFile( inpFile=inpFile, returnType="structured" )

    # ------------------------------------------------- #
    # --- [2] store in grid                         --- #
    # ------------------------------------------------- #

    import nkBasicAlgs.store__inGrid3D as sig
    BField_   = sig.store__inGrid3D( Data=BField, x1MinMaxNum=const["x1MinMaxNum"], \
                                     x2MinMaxNum=const["x2MinMaxNum"], \
                                     x3MinMaxNum=const["x3MinMaxNum"],  )
    source[:,:,:,3:6] = BField_[:,:,:,3:6]

    # ------------------------------------------------- #
    # --- [3] save as source.dat                    --- #
    # ------------------------------------------------- #

    import nkUtilities.save__pointFile as spf
    outFile   = "dat/source.dat"
    spf.save__pointFile( outFile=outFile, Data=source )


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    load__field()

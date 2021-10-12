# read_input_script.py:
# Reads the names of the emission line provided (yes or no) in the input.txt file. 
# Reads the fits 2D flux and error images of emission lines and save them in 2D arrays. 
# These new 2D arrays are bigger that the first ones because the center of the nebula 
# is moved to the center of the new map and some rows/columns of pixels have been added.
# (C) Stavros Akras

import numpy as np
from astropy.io.fits import getdata
from satellite import fluxes as flx



def read_input_images(line_name,irequested,param,flux_or_error):
    
##########################################################################
### extra rows/columns of pixels are create and then added in the raw maps
### the number of rows and columns are given by the user.
##########################################################################
    z1=np.zeros((param[3],param[7]))  # above
    z2=np.zeros((param[4],param[7]))  # below
    z3=np.zeros((param[8],param[5]))  # left 
    z4=np.zeros((param[8],param[6]))  #right
    
    flux = flx.Flux2D()
    
    if flux_or_error=="fluxes":
        for linet in [l for l in zip(line_name, irequested) if l[0].endswith('s')]:
            linem, requested = linet[0], linet[1]
            if linem == 'HI_6563s' and requested=='yes':
                dataHa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHa=np.nan_to_num(dataHa)
                addupdown1=np.concatenate((z2,dataHa,z1), axis=0)
                flux.Ha_6563=np.concatenate((z3,addupdown1,z4), axis=1)
                sx=hdr['NAXIS1']
                sy=hdr['NAXIS2']
            elif linem == 'HI_4861s' and requested=='yes':
                dataHb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHb=np.nan_to_num(dataHb)
                addupdown1=np.concatenate((z2,dataHb,z1), axis=0)
                flux.Hb_4861=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HI_4861s' and requested=='no':
                dataHb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHb,z1), axis=0)
                flux.Hb_4861=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HI_4340s' and requested=='yes':
                dataHg, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHg=np.nan_to_num(dataHg)
                addupdown1=np.concatenate((z2,dataHg,z1), axis=0)
                flux.Hg_4340=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HI_4340s' and requested=='no':
                dataHg=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHg,z1), axis=0)
                flux.Hg_4340=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HI_4101s' and requested=='yes':
                dataHd, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHd=np.nan_to_num(dataHd)
                addupdown1=np.concatenate((z2,dataHd,z1), axis=0)
                flux.Hd_4101=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HI_4101s' and requested=='no':
                dataHd=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHd,z1), axis=0)
                flux.Hd_4101=np.concatenate((z3,addupdown1,z4), axis=1)

            elif linem == 'HeI_5876s' and requested=='yes':
                dataHeIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHeIa=np.nan_to_num(dataHeIa)
                addupdown1=np.concatenate((z2,dataHeIa,z1), axis=0)
                flux.HeIa_5876=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HeI_5876s' and requested=='no':
                dataHeIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIa,z1), axis=0)
                flux.HeIa_5876=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HeI_6678s' and requested=='yes':
                dataHeIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHeIb=np.nan_to_num(dataHeIb)
                addupdown1=np.concatenate((z2,dataHeIb,z1), axis=0)
                flux.HeIb_6678=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HeI_6678s' and requested=='no':
                dataHeIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIb,z1), axis=0)
                flux.HeIb_6678=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HeII_4686s' and requested=='yes':    
                dataHeIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHeIIa=np.nan_to_num(dataHeIIa)
                addupdown1=np.concatenate((z2,dataHeIIa,z1), axis=0)
                flux.HeIIa_4686=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HeII_4686s' and requested=='no':
                dataHeIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIIa,z1), axis=0)
                flux.HeIIa_4686=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HeII_5412s' and requested=='yes':
                dataHeIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataHeIIb=np.nan_to_num(dataHeIIb)
                addupdown1=np.concatenate((z2,dataHeIIb,z1), axis=0)
                flux.HeIIb_5412=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'HeII_5412s' and requested=='no':
                dataHeIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIIb,z1), axis=0)
                flux.HeIIb_5412=np.concatenate((z3,addupdown1,z4), axis=1)
            
            elif linem == 'N2_5755s' and requested=='yes':
                dataNIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataNIIa=np.nan_to_num(dataNIIa)
                addupdown1=np.concatenate((z2,dataNIIa,z1), axis=0)
                flux.NIIa_5755=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'N2_5755s' and requested=='no':
                dataNIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIa,z1), axis=0)
                flux.NIIa_5755=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'N2_6548s' and requested=='yes':
                dataNIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataNIIb=np.nan_to_num(dataNIIb)
                addupdown1=np.concatenate((z2,dataNIIb,z1), axis=0)
                flux.NIIb_6548=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'N2_6548s' and requested=='no':
                dataNIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIb,z1), axis=0)
                flux.NIIb_6548=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'N2_6583s' and requested=='yes':
                dataNIIc, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataNIIc=np.nan_to_num(dataNIIc)
                addupdown1=np.concatenate((z2,dataNIIc,z1), axis=0)
                flux.NIIc_6584=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'N2_6583s' and requested=='no':
                dataNIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIc,z1), axis=0)
                flux.NIIc_6584=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'N1_5199s' and requested=='yes':
                dataNI, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataNI=np.nan_to_num(dataNI)
                addupdown1=np.concatenate((z2,dataNI,z1), axis=0)
                flux.NI_5199=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'N1_5199s' and requested=='no':
                dataNI=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNI,z1), axis=0)
                flux.NI_5199=np.concatenate((z3,addupdown1,z4), axis=1)
            
            elif linem == 'O3_4363s' and requested=='yes':
                dataOIIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIIIa=np.nan_to_num(dataOIIIa)
                addupdown1=np.concatenate((z2,dataOIIIa,z1), axis=0)
                flux.OIIIa_4363=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O3_4363s' and requested=='no':
                dataOIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIa,z1), axis=0)
                flux.OIIIa_4363=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O3_4959s' and requested=='yes':
                dataOIIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIIIb=np.nan_to_num(dataOIIIb)
                addupdown1=np.concatenate((z2,dataOIIIb,z1), axis=0)
                flux.OIIIb_4959=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O3_4959s' and requested=='no':
                dataOIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIb,z1), axis=0)
                flux.OIIIb_4959=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O3_5007s' and requested=='yes':
                dataOIIIc, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIIIc=np.nan_to_num(dataOIIIc)
                addupdown1=np.concatenate((z2,dataOIIIc,z1), axis=0)
                flux.OIIIc_5007=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O3_5007s' and requested=='no':
                dataOIIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIc,z1), axis=0)
                flux.OIIIc_5007=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_3727s' and requested=='yes':
                dataOIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIIa=np.nan_to_num(dataOIIa)
                addupdown1=np.concatenate((z2,dataOIIa,z1), axis=0)
                flux.OIIa_3727=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_3727s' and requested=='no':
                dataOIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIa,z1), axis=0)
                flux.OIIa_3727=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_3729s' and requested=='yes':
                dataOIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIIb=np.nan_to_num(dataOIIb)
                addupdown1=np.concatenate((z2,dataOIIb,z1), axis=0)
                flux.OIIb_3729=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_3729s' and requested=='no':
                dataOIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIb,z1), axis=0)
                flux.OIIb_3729=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_7320s' and requested=='yes':
                dataOIIc, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIIc=np.nan_to_num(dataOIIc)
                addupdown1=np.concatenate((z2,dataOIIc,z1), axis=0)
                flux.OIIc_7320=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_7320s' and requested=='no':
                dataOIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIc,z1), axis=0)
                flux.OIIc_7320=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_7330s' and requested=='yes':
                dataOIId, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIId=np.nan_to_num(dataOIId)
                addupdown1=np.concatenate((z2,dataOIId,z1), axis=0)
                flux.OIId_7330=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O2_7330s' and requested=='no':
                dataOIId=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIId,z1), axis=0)
                flux.OIId_7330=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O1_5577s' and requested=='yes':
                dataOIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIa=np.nan_to_num(dataOIa)
                addupdown1=np.concatenate((z2,dataOIa,z1), axis=0)
                flux.OIa_5577=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O1_5577s' and requested=='no':
                dataOIa=np.zeros((sy, sx))   
                addupdown1=np.concatenate((z2,dataOIa,z1), axis=0)
                flux.OIa_5577=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O1_6300s' and requested=='yes':
                dataOIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIb=np.nan_to_num(dataOIb)
                addupdown1=np.concatenate((z2,dataOIb,z1), axis=0)
                flux.OIb_6300=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O1_6300s' and requested=='no':
                dataOIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIb,z1), axis=0)
                flux.OIb_6300=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O1_6363s' and requested=='yes':
                dataOIc, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataOIc=np.nan_to_num(dataOIc)
                addupdown1=np.concatenate((z2,dataOIc,z1), axis=0)
                flux.OIc_6363=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'O1_6363s' and requested=='no':
                dataOIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIc,z1), axis=0)
                flux.OIc_6363=np.concatenate((z3,addupdown1,z4), axis=1)
                    
            elif linem == 'S2_6716s' and requested=='yes':
                dataSIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataSIIa=np.nan_to_num(dataSIIa)
                addupdown1=np.concatenate((z2,dataSIIa,z1), axis=0)
                flux.SIIa_6716=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'S2_6716s' and requested=='no':
                dataSIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIa,z1), axis=0)
                flux.SIIa_6716=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'S2_6731s' and requested=='yes':
                dataSIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataSIIb=np.nan_to_num(dataSIIb)
                addupdown1=np.concatenate((z2,dataSIIb,z1), axis=0)
                flux.SIIb_6731=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'S2_6731s' and requested=='no':
                dataSIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIb,z1), axis=0)
                flux.SIIb_6731=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'S3_6312s' and requested=='yes':
                dataSIIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataSIIIa=np.nan_to_num(dataSIIIa)
                addupdown1=np.concatenate((z2,dataSIIIa,z1), axis=0)
                flux.SIIIa_6312=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'S3_6312s' and requested=='no':
                dataSIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIIa,z1), axis=0)
                flux.SIIIa_6312=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'S3_9069s' and requested=='yes':
                dataSIIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataSIIIb=np.nan_to_num(dataSIIIb)
                addupdown1=np.concatenate((z2,dataSIIIb,z1), axis=0)
                flux.SIIIb_9069=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'S3_9069s' and requested=='no':
                dataSIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIIb,z1), axis=0)
                flux.SIIIb_9069=np.concatenate((z3,addupdown1,z4), axis=1)
                
            elif linem == 'Cl3_5517s' and requested=='yes':
                dataClIIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataClIIIa=np.nan_to_num(dataClIIIa)
                addupdown1=np.concatenate((z2,dataClIIIa,z1), axis=0)
                flux.ClIIIa_5517=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Cl3_5517s' and requested=='no':
                dataClIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataClIIIa,z1), axis=0)
                flux.ClIIIa_5517=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Cl3_5538s' and requested=='yes':
                dataClIIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataClIIIb=np.nan_to_num(dataClIIIb)
                addupdown1=np.concatenate((z2,dataClIIIb,z1), axis=0)
                flux.ClIIIb_5538=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Cl3_5538s' and requested=='no':
                dataClIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataClIIIb,z1), axis=0)
                flux.ClIIIb_5538=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ar3_7136s' and requested=='yes':
                dataArIII, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataArIII=np.nan_to_num(dataArIII)
                addupdown1=np.concatenate((z2,dataArIII,z1), axis=0)
                flux.ArIII_7136=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ar3_7136s' and requested=='no':
                dataArIII=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIII,z1), axis=0)
                flux.ArIII_7136=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ar4_4712s' and requested=='yes':
                dataArIVa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataArIVa=np.nan_to_num(dataArIVa)
                addupdown1=np.concatenate((z2,dataArIVa,z1), axis=0)
                flux.ArIVa_4712=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ar4_4712s' and requested=='no':
                dataArIVa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIVa,z1), axis=0)
                flux.ArIVa_4712=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ar4_4740s' and requested=='yes':
                dataArIVb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataArIVb=np.nan_to_num(dataArIVb)
                addupdown1=np.concatenate((z2,dataArIVb,z1), axis=0)
                flux.ArIVb_4740=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ar4_4740s' and requested=='no':
                dataArIVb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIVb,z1), axis=0)
                flux.ArIVb_4740=np.concatenate((z3,addupdown1,z4), axis=1)

            elif linem == 'C1_8727s' and requested=='yes':
                dataCI, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataCI=np.nan_to_num(dataCI)
                addupdown1=np.concatenate((z2,dataCI,z1), axis=0)
                flux.CI_8727=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'C1_8727s' and requested=='no':
                dataCI=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataCI,z1), axis=0)
                flux.CI_8727=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'C2_6461s' and requested=='yes':
                dataCII, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataCII=np.nan_to_num(dataCII)
                addupdown1=np.concatenate((z2,dataCII,z1), axis=0)
                flux.CII_6461=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'C2_6461s' and requested=='no':
                dataCII=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataCII,z1), axis=0)
                flux.CII_6461=np.concatenate((z3,addupdown1,z4), axis=1)
    
            elif linem == 'Ne3_3868s' and requested=='yes':
                dataNeIIIa, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataNeIIIa=np.nan_to_num(dataNeIIIa)
                addupdown1=np.concatenate((z2,dataNeIIIa,z1), axis=0)
                flux.NeIIIa_3868=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ne3_3868s' and requested=='no':
                dataNeIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNeIIIa,z1), axis=0)
                flux.NeIIIa_3868=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ne3_3967s' and requested=='yes':
                dataNeIIIb, hdr = getdata('image_data/'+linem+'.fit', 0, header=True)
                dataNeIIIb=np.nan_to_num(dataNeIIIb)
                addupdown1=np.concatenate((z2,dataNeIIIb,z1), axis=0)
                flux.NeIIIb_3967=np.concatenate((z3,addupdown1,z4), axis=1)
            elif linem == 'Ne3_3967s' and requested=='no':
                dataNeIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNeIIIb,z1), axis=0)
                flux.NeIIIb_3967=np.concatenate((z3,addupdown1,z4), axis=1)
        
        return flux,hdr

    if flux_or_error=="errors":
        for linet in [l for l in zip(line_name, irequested) if l[0].endswith('e')]:
            linem, requested = linet[0], linet[1]
            if linem == 'HI_6563e' and requested == 'yes':
                dataHa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHa = np.nan_to_num(dataHa)
                addupdown1 = np.concatenate((z2, dataHa, z1), axis=0)
                flux.Ha_6563 = np.concatenate((z3, addupdown1, z4), axis=1)
                sx = hdr['NAXIS1']
                sy = hdr['NAXIS2']
            elif linem == 'HI_4861e' and requested == 'yes':
                dataHb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHb = np.nan_to_num(dataHb)
                addupdown1 = np.concatenate((z2, dataHb, z1), axis=0)
                flux.Hb_4861 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HI_4861e' and requested == 'no':
                dataHb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataHb, z1), axis=0)
                flux.Hb_4861 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HI_4340e' and requested == 'yes':
                dataHg, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHg = np.nan_to_num(dataHg)
                addupdown1 = np.concatenate((z2, dataHg, z1), axis=0)
                flux.Hg_4340 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HI_4340e' and requested == 'no':
                dataHg = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataHg, z1), axis=0)
                flux.Hg_4340 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HI_4101e' and requested == 'yes':
                dataHd, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHd = np.nan_to_num(dataHd)
                addupdown1 = np.concatenate((z2, dataHd, z1), axis=0)
                flux.Hd_4101 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HI_4101e' and requested == 'no':
                dataHd = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataHd, z1), axis=0)
                flux.Hd_4101 = np.concatenate((z3, addupdown1, z4), axis=1)

            elif linem == 'HeI_5876e' and requested == 'yes':
                dataHeIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHeIa = np.nan_to_num(dataHeIa)
                addupdown1 = np.concatenate((z2, dataHeIa, z1), axis=0)
                flux.HeIa_5876 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HeI_5876e' and requested == 'no':
                dataHeIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataHeIa, z1), axis=0)
                flux.HeIa_5876 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HeI_6678e' and requested == 'yes':
                dataHeIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHeIb = np.nan_to_num(dataHeIb)
                addupdown1 = np.concatenate((z2, dataHeIb, z1), axis=0)
                flux.HeIb_6678 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HeI_6678e' and requested == 'no':
                dataHeIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataHeIb, z1), axis=0)
                flux.HeIb_6678 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HeII_4686e' and requested == 'yes':
                dataHeIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHeIIa = np.nan_to_num(dataHeIIa)
                addupdown1 = np.concatenate((z2, dataHeIIa, z1), axis=0)
                flux.HeIIa_4686 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HeII_4686e' and requested == 'no':
                dataHeIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataHeIIa, z1), axis=0)
                flux.HeIIa_4686 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HeII_5412e' and requested == 'yes':
                dataHeIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataHeIIb = np.nan_to_num(dataHeIIb)
                addupdown1 = np.concatenate((z2, dataHeIIb, z1), axis=0)
                flux.HeIIb_5412 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'HeII_5412e' and requested == 'no':
                dataHeIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataHeIIb, z1), axis=0)
                flux.HeIIb_5412 = np.concatenate((z3, addupdown1, z4), axis=1)

            elif linem == 'N2_5755e' and requested == 'yes':
                dataNIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataNIIa = np.nan_to_num(dataNIIa)
                addupdown1 = np.concatenate((z2, dataNIIa, z1), axis=0)
                flux.NIIa_5755 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'N2_5755e' and requested == 'no':
                dataNIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataNIIa, z1), axis=0)
                flux.NIIa_5755 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'N2_6548e' and requested == 'yes':
                dataNIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataNIIb = np.nan_to_num(dataNIIb)
                addupdown1 = np.concatenate((z2, dataNIIb, z1), axis=0)
                flux.NIIb_6548 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'N2_6548e' and requested == 'no':
                dataNIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataNIIb, z1), axis=0)
                flux.NIIb_6548 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'N2_6583e' and requested == 'yes':
                dataNIIc, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataNIIc = np.nan_to_num(dataNIIc)
                addupdown1 = np.concatenate((z2, dataNIIc, z1), axis=0)
                flux.NIIc_6584 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'N2_6583e' and requested == 'no':
                dataNIIc = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataNIIc, z1), axis=0)
                flux.NIIc_6584 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'N1_5199e' and requested == 'yes':
                dataNI, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataNI = np.nan_to_num(dataNI)
                addupdown1 = np.concatenate((z2, dataNI, z1), axis=0)
                flux.NI_5199 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'N1_5199e' and requested == 'no':
                dataNI = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataNI, z1), axis=0)
                flux.NI_5199 = np.concatenate((z3, addupdown1, z4), axis=1)

            elif linem == 'O3_4363e' and requested == 'yes':
                dataOIIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIIIa = np.nan_to_num(dataOIIIa)
                addupdown1 = np.concatenate((z2, dataOIIIa, z1), axis=0)
                flux.OIIIa_4363 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O3_4363e' and requested == 'no':
                dataOIIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIIIa, z1), axis=0)
                flux.OIIIa_4363 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O3_4959e' and requested == 'yes':
                dataOIIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIIIb = np.nan_to_num(dataOIIIb)
                addupdown1 = np.concatenate((z2, dataOIIIb, z1), axis=0)
                flux.OIIIb_4959 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O3_4959e' and requested == 'no':
                dataOIIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIIIb, z1), axis=0)
                flux.OIIIb_4959 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O3_5007e' and requested == 'yes':
                dataOIIIc, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIIIc = np.nan_to_num(dataOIIIc)
                addupdown1 = np.concatenate((z2, dataOIIIc, z1), axis=0)
                flux.OIIIc_5007 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O3_5007e' and requested == 'no':
                dataOIIIc = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIIIc, z1), axis=0)
                flux.OIIIc_5007 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_3727e' and requested == 'yes':
                dataOIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIIa = np.nan_to_num(dataOIIa)
                addupdown1 = np.concatenate((z2, dataOIIa, z1), axis=0)
                flux.OIIa_3727 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_3727e' and requested == 'no':
                dataOIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIIa, z1), axis=0)
                flux.OIIa_3727 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_3729e' and requested == 'yes':
                dataOIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIIb = np.nan_to_num(dataOIIb)
                addupdown1 = np.concatenate((z2, dataOIIb, z1), axis=0)
                flux.OIIb_3729 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_3729e' and requested == 'no':
                dataOIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIIb, z1), axis=0)
                flux.OIIb_3729 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_7320e' and requested == 'yes':
                dataOIIc, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIIc = np.nan_to_num(dataOIIc)
                addupdown1 = np.concatenate((z2, dataOIIc, z1), axis=0)
                flux.OIIc_7320 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_7320e' and requested == 'no':
                dataOIIc = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIIc, z1), axis=0)
                flux.OIIc_7320 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_7330e' and requested == 'yes':
                dataOIId, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIId = np.nan_to_num(dataOIId)
                addupdown1 = np.concatenate((z2, dataOIId, z1), axis=0)
                flux.OIId_7330 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O2_7330e' and requested == 'no':
                dataOIId = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIId, z1), axis=0)
                flux.OIId_7330 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O1_5577e' and requested == 'yes':
                dataOIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIa = np.nan_to_num(dataOIa)
                addupdown1 = np.concatenate((z2, dataOIa, z1), axis=0)
                flux.OIa_5577 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O1_5577e' and requested == 'no':
                dataOIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIa, z1), axis=0)
                flux.OIa_5577 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O1_6300e' and requested == 'yes':
                dataOIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIb = np.nan_to_num(dataOIb)
                addupdown1 = np.concatenate((z2, dataOIb, z1), axis=0)
                flux.OIb_6300 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O1_6300e' and requested == 'no':
                dataOIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIb, z1), axis=0)
                flux.OIb_6300 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O1_6363e' and requested == 'yes':
                dataOIc, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataOIc = np.nan_to_num(dataOIc)
                addupdown1 = np.concatenate((z2, dataOIc, z1), axis=0)
                flux.OIc_6363 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'O1_6363e' and requested == 'no':
                dataOIc = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataOIc, z1), axis=0)
                flux.OIc_6363 = np.concatenate((z3, addupdown1, z4), axis=1)

            elif linem == 'S2_6716e' and requested == 'yes':
                dataSIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataSIIa = np.nan_to_num(dataSIIa)
                addupdown1 = np.concatenate((z2, dataSIIa, z1), axis=0)
                flux.SIIa_6716 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'S2_6716e' and requested == 'no':
                dataSIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataSIIa, z1), axis=0)
                flux.SIIa_6716 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'S2_6731e' and requested == 'yes':
                dataSIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataSIIb = np.nan_to_num(dataSIIb)
                addupdown1 = np.concatenate((z2, dataSIIb, z1), axis=0)
                flux.SIIb_6731 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'S2_6731e' and requested == 'no':
                dataSIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataSIIb, z1), axis=0)
                flux.SIIb_6731 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'S3_6312e' and requested == 'yes':
                dataSIIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataSIIIa = np.nan_to_num(dataSIIIa)
                addupdown1 = np.concatenate((z2, dataSIIIa, z1), axis=0)
                flux.SIIIa_6312 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'S3_6312e' and requested == 'no':
                dataSIIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataSIIIa, z1), axis=0)
                flux.SIIIa_6312 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'S3_9069e' and requested == 'yes':
                dataSIIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataSIIIb = np.nan_to_num(dataSIIIb)
                addupdown1 = np.concatenate((z2, dataSIIIb, z1), axis=0)
                flux.SIIIb_9069 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'S3_9069e' and requested == 'no':
                dataSIIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataSIIIb, z1), axis=0)
                flux.SIIIb_9069 = np.concatenate((z3, addupdown1, z4), axis=1)

            elif linem == 'Cl3_5517e' and requested == 'yes':
                dataClIIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataClIIIa = np.nan_to_num(dataClIIIa)
                addupdown1 = np.concatenate((z2, dataClIIIa, z1), axis=0)
                flux.ClIIIa_5517 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Cl3_5517e' and requested == 'no':
                dataClIIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataClIIIa, z1), axis=0)
                flux.ClIIIa_5517 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Cl3_5538e' and requested == 'yes':
                dataClIIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataClIIIb = np.nan_to_num(dataClIIIb)
                addupdown1 = np.concatenate((z2, dataClIIIb, z1), axis=0)
                flux.ClIIIb_5538 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Cl3_5538e' and requested == 'no':
                dataClIIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataClIIIb, z1), axis=0)
                flux.ClIIIb_5538 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ar3_7136e' and requested == 'yes':
                dataArIII, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataArIII = np.nan_to_num(dataArIII)
                addupdown1 = np.concatenate((z2, dataArIII, z1), axis=0)
                flux.ArIII_7136 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ar3_7136e' and requested == 'no':
                dataArIII = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataArIII, z1), axis=0)
                flux.ArIII_7136 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ar4_4712e' and requested == 'yes':
                dataArIVa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataArIVa = np.nan_to_num(dataArIVa)
                addupdown1 = np.concatenate((z2, dataArIVa, z1), axis=0)
                flux.ArIVa_4712 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ar4_4712e' and requested == 'no':
                dataArIVa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataArIVa, z1), axis=0)
                flux.ArIVa_4712 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ar4_4740e' and requested == 'yes':
                dataArIVb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataArIVb = np.nan_to_num(dataArIVb)
                addupdown1 = np.concatenate((z2, dataArIVb, z1), axis=0)
                flux.ArIVb_4740 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ar4_4740e' and requested == 'no':
                dataArIVb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataArIVb, z1), axis=0)
                flux.ArIVb_4740 = np.concatenate((z3, addupdown1, z4), axis=1)

            elif linem == 'C1_8727e' and requested == 'yes':
                dataCI, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataCI = np.nan_to_num(dataCI)
                addupdown1 = np.concatenate((z2, dataCI, z1), axis=0)
                flux.CI_8727 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'C1_8727e' and requested == 'no':
                dataCI = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataCI, z1), axis=0)
                flux.CI_8727 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'C2_6461e' and requested == 'yes':
                dataCII, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataCII = np.nan_to_num(dataCII)
                addupdown1 = np.concatenate((z2, dataCII, z1), axis=0)
                flux.CII_6461 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'C2_6461e' and requested == 'no':
                dataCII = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataCII, z1), axis=0)
                flux.CII_6461 = np.concatenate((z3, addupdown1, z4), axis=1)

            elif linem == 'Ne3_3868e' and requested == 'yes':
                dataNeIIIa, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataNeIIIa = np.nan_to_num(dataNeIIIa)
                addupdown1 = np.concatenate((z2, dataNeIIIa, z1), axis=0)
                flux.NeIIIa_3868 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ne3_3868e' and requested == 'no':
                dataNeIIIa = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataNeIIIa, z1), axis=0)
                flux.NeIIIa_3868 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ne3_3967e' and requested == 'yes':
                dataNeIIIb, hdr = getdata(
                    'image_data/'+linem+'.fit', 0, header=True)
                dataNeIIIb = np.nan_to_num(dataNeIIIb)
                addupdown1 = np.concatenate((z2, dataNeIIIb, z1), axis=0)
                flux.NeIIIb_3967 = np.concatenate((z3, addupdown1, z4), axis=1)
            elif linem == 'Ne3_3967e' and requested == 'no':
                dataNeIIIb = np.zeros((sy, sx))
                addupdown1 = np.concatenate((z2, dataNeIIIb, z1), axis=0)
                flux.NeIIIb_3967 = np.concatenate((z3, addupdown1, z4), axis=1)

        return flux,hdr
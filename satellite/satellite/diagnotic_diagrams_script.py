# diagnotic_diagrams_script.py:
#
# (C) Stavros Akras

from __future__ import print_function

import numpy as np

from satellite import line_ratio_script as lrs
from satellite import positive_values_script as pvs
from satellite import mean_value_script as mvs
from satellite import std_value_script as svs

def generic_plot(numerator1, denominator1, numerator2, denominator2, sizex, sizey, lineratio_index, label1, label2):

    with open('general_output_file.txt', 'a') as f10:
        data_ratio = lrs.lineratio(numerator1, denominator1, sizex, sizey, lineratio_index)
        #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
        data_pos_ratio1 = pvs.positivevalues(data_ratio, numerator1, denominator1,
                                             sizex, sizey, 0)
        mean_ratio = mvs.meanvalue(data_pos_ratio1, -10)
        std_ratio = svs.stdvalue(data_pos_ratio1, mean_ratio, -10)
        print("log({:}) {:.3f} sigma {:.3f}".format(
            label1, mean_ratio, std_ratio), file=f10)

        data_ratio = lrs.lineratio(numerator2, denominator2, sizex, sizey, lineratio_index)
        #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
        data_pos_ratio2 = pvs.positivevalues(data_ratio, numerator2, denominator2,
                                             sizex, sizey, 0)
        mean_ratio = mvs.meanvalue(data_pos_ratio2, -10)
        std_ratio = svs.stdvalue(data_pos_ratio2, mean_ratio, -10)
        print("log({:}) {:.3f} sigma {:.3f}".format(
            label2, mean_ratio, std_ratio), file=f10)

    maxper = 75
    minper = 25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

    #    per25n=np.percentile(data_pos_ratio1,minper)
    #    per75n=np.percentile(data_pos_ratio1,maxper)
    #    per25s=np.percentile(data_pos_ratio2,minper)
    #    per75s=np.percentile(data_pos_ratio2,maxper)

    datax = [-99] * (sizex * sizey)
    datay = [-99] * (sizex * sizey)
    k = 0
    for j in range(sizex):
        for i in range(sizey):

            datax.append(-99)
            datay.append(-99)

            if (numerator1[i, j] == 0 or denominator1[i, j] == 0):
                datay[k] = -99

            if denominator1[i, j] < 0:
                datay[k] = -99

            if numerator1[i, j] > 0 and denominator1[i, j] > 0:
                datay[k] = data_ratio1[i, j]

            if (numerator2[i, j] == 0 or denominator2[i, j] == 0):
                datax[k] = -99

            if denominator2[i, j] < 0:
                datax[k] = -99

            if numerator2[i, j] > 0 and denominator2[i, j] > 0:
                datax[k] = data_ratio2[i, j]

            k = k + 1

    xpos = [-99] * len(datax)
    ypos = [-99] * len(datax)
    l = 0
    for j in range(0, len(datax)):
        if (datay[j] > -98 and datax[j] > -98):
            xpos.append(-99)
            ypos.append(-99)
            xpos[l] = datax[j]
            ypos[l] = datay[j]
            l = l + 1


#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT

    return xpos, ypos


def NIISIIplot(numerator1, denominator1, numerator2, denominator2, sizex,
               sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 1, sizex, sizey, "Ha/[N II] 6548+6584", "Ha/[S II] 6716+6731")

def OIIINIIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[N II] 6584/Ha")

def OIIISIIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[S II] 6716+6731/Ha")

def OIIIHeII4686plot(numerator1, denominator1, numerator2, denominator2, sizex,
                     sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "He II 4686/Hb")

def OIIIHeII5412plot(numerator1, denominator1, numerator2, denominator2, sizex,
                     sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "He II 5412/Hb")

def OIIINIplot(numerator1, denominator1, numerator2, denominator2, sizex,
               sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[N I] 5199/Hb")

def OIIIOIplot(numerator1, denominator1, numerator2, denominator2, sizex,
               sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[O I] 6300/Ha")

def OIIINII5755plot(numerator1, denominator1, numerator2, denominator2, sizex,
                    sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[N II] 5755/Ha")

def SII67166731HaSIIplot(numerator1, denominator1, numerator2, denominator2,
                         sizex, sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 2, sizex, sizey, "[S II] 6716/6731", "[Ha/[S II] 6716+6731")

def SII67166731HaNIIplot(numerator1, denominator1, numerator2, denominator2,
                         sizex, sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 2, sizex, sizey, "[S II] 6716/6731", "[Ha/[N II] 6548+6584")

def NIISIISIIOIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                    sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[N II] 6548+6584/[S II] 6716+6731", "[S II] 6716+6731/[O I] 6300+6363")

def NIISIINIIOIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                    sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[N II] 6548+6584/[S II] 6716+6731", "[N II] 658+6584/[O I] 6300+6363")

def OIIINII6584OIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                      sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[N II] 6584/[O I] 6300")

def OIIIOII732030OIIIplot(numerator1, denominator1, numerator2, denominator2,
                          sizex, sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[O II] 7320+7330/[O III] 5007")

def OIIIOII372729OIIIplot(numerator1, denominator1, numerator2, denominator2,
                          sizex, sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[O II] 3727+3729/[O III] 5007")

def OIIIOII372729Hbplot(numerator1, denominator1, numerator2, denominator2,
                        sizex, sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[O II] 3727+3729/Hb")

def OIIIOII732030Haplot(numerator1, denominator1, numerator2, denominator2,
                        sizex, sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[O II] 7320+7330/Ha")

def OIIIHeIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "He I 5876/Ha")

def OIIIArIIIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                  sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[Ar III] 7136/Ha")

def OIIIOIplot(numerator1, denominator1, numerator2, denominator2, sizex,
               sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[O III] 5007/[O I] 6300")

def OIIININIIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                  sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[N I] 5199/[N II] 6584")

def OIIINeIIIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                  sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[Ne III] 3869/Hb")

def OIIIOIIIHgplot(numerator1, denominator1, numerator2, denominator2, sizex,
                   sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[O III] 4363/Hg")

def OIIIArIVplot(numerator1, denominator1, numerator2, denominator2, sizex,
                 sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[Ar IV] 4712/4740")

def OIIICIplot(numerator1, denominator1, numerator2, denominator2, sizex,
               sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[C I] 8727/Ha")

def OIIICIIplot(numerator1, denominator1, numerator2, denominator2, sizex,
                sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O III] 5007/Hb", "[C II] 6461/Ha")

def OICIplot(numerator1, denominator1, numerator2, denominator2, sizex, sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[O I] 6300/Ha", "[C I] 8727/Ha")

def ArIVHeII4686plot(numerator1, denominator1, numerator2, denominator2, sizex,
                     sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[Ar IV] 4712+4740/Hb", "He II 4686/Hb")

def ArIVHeII5412plot(numerator1, denominator1, numerator2, denominator2, sizex,
                     sizey):
    return generic_plot(numerator1, denominator1, numerator2, denominator2, 0, sizex, sizey, "[Ar IV] 4712+4740/Hb", "He II 5412/Hb")
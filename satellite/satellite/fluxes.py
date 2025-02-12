import numpy as np


class Flux2D:

    def __init__(self, sizex=0, sizey=0):
        if sizex or sizey:
            self.Ha_6563 = np.zeros((sizex, sizey))
            self.Hb_4861 = np.zeros((sizex, sizey))
            self.Hg_4340 = np.zeros((sizex, sizey))
            self.Hd_4101 = np.zeros((sizex, sizey))
            self.HeIa_5876 = np.zeros((sizex, sizey))
            self.HeIb_6678 = np.zeros((sizex, sizey))
            self.HeIIa_4686 = np.zeros((sizex, sizey))
            self.HeIIb_5412 = np.zeros((sizex, sizey))
            self.NIIa_5755 = np.zeros((sizex, sizey))
            self.NIIb_6548 = np.zeros((sizex, sizey))
            self.NIIc_6584 = np.zeros((sizex, sizey))
            self.NI_5199 = np.zeros((sizex, sizey))
            self.OIIIa_4363 = np.zeros((sizex, sizey))
            self.OIIIb_4959 = np.zeros((sizex, sizey))
            self.OIIIc_5007 = np.zeros((sizex, sizey))
            self.OIIa_3727 = np.zeros((sizex, sizey))
            self.OIIb_3729 = np.zeros((sizex, sizey))
            self.OIIc_7320 = np.zeros((sizex, sizey))
            self.OIId_7330 = np.zeros((sizex, sizey))
            self.OIa_5577 = np.zeros((sizex, sizey))
            self.OIb_6300 = np.zeros((sizex, sizey))
            self.OIc_6363 = np.zeros((sizex, sizey))
            self.SIIa_6716 = np.zeros((sizex, sizey))
            self.SIIb_6731 = np.zeros((sizex, sizey))
            self.SIIIa_6312 = np.zeros((sizex, sizey))
            self.SIIIb_9069 = np.zeros((sizex, sizey))
            self.ClIIIa_5517 = np.zeros((sizex, sizey))
            self.ClIIIb_5538 = np.zeros((sizex, sizey))
            self.ArIII_7136 = np.zeros((sizex, sizey))
            self.ArIVa_4712 = np.zeros((sizex, sizey))
            self.ArIVb_4740 = np.zeros((sizex, sizey))
            self.CI_8727 = np.zeros((sizex, sizey))
            self.CII_6461 = np.zeros((sizex, sizey))
            self.NeIIIa_3868 = np.zeros((sizex, sizey))
            self.NeIIIb_3967 = np.zeros((sizex, sizey))
        else:
            self.Ha_6563 = []
            self.Hb_4861 = []
            self.Hg_4340 = []
            self.Hd_4101 = []
            self.HeIa_5876 = []
            self.HeIb_6678 = []
            self.HeIIa_4686 = []
            self.HeIIb_5412 = []
            self.NIIa_5755 = []
            self.NIIb_6548 = []
            self.NIIc_6584 = []
            self.NI_5199 = []
            self.OIIIa_4363 = []
            self.OIIIb_4959 = []
            self.OIIIc_5007 = []
            self.OIIa_3727 = []
            self.OIIb_3729 = []
            self.OIIc_7320 = []
            self.OIId_7330 = []
            self.OIa_5577 = []
            self.OIb_6300 = []
            self.OIc_6363 = []
            self.SIIa_6716 = []
            self.SIIb_6731 = []
            self.SIIIa_6312 = []
            self.SIIIb_9069 = []
            self.ClIIIa_5517 = []
            self.ClIIIb_5538 = []
            self.ArIII_7136 = []
            self.ArIVa_4712 = []
            self.ArIVb_4740 = []
            self.CI_8727 = []
            self.CII_6461 = []
            self.NeIIIa_3868 = []
            self.NeIIIb_3967 = []


class TeNe:

    def __init__(self, sizex=0, sizey=0):
        if sizex !=0 or sizey !=0:
            if sizex != 0 and sizey == 0: sizey = sizex
            self.NIISII = np.zeros((sizex, sizey))
            self.OISII = np.zeros((sizex, sizey))
            self.OIISII = np.zeros((sizex, sizey))
            self.OIIISII = np.zeros((sizex, sizey))
            self.SIIISII = np.zeros((sizex, sizey))
            self.OIIOII = np.zeros((sizex, sizey))
            self.NIIOII = np.zeros((sizex, sizey))
            self.OIOII = np.zeros((sizex, sizey))
            self.OIIIClIII = np.zeros((sizex, sizey))
            self.SIIIClIII = np.zeros((sizex, sizey))
            self.OIIIArVI = np.zeros((sizex, sizey))
            self.SIIIArVI = np.zeros((sizex, sizey))
            self.NIIClIII = np.zeros((sizex, sizey))
        else:
            self.NIISII = []
            self.OISII = []
            self.OIISII = []
            self.OIIISII = []
            self.SIIISII = []
            self.OIIOII = []
            self.NIIOII = []
            self.OIOII = []
            self.OIIIClIII = []
            self.SIIIClIII = []
            self.OIIIArVI = []
            self.SIIIArVI = []
            self.NIIClIII = []
    
    def append(self, value, times):
        lst = [value] * times
        for atr in [self.NIISII, self.OISII, self.OIISII, self.OIIISII, self.SIIISII, self.OIIOII, self.NIIOII, self.OIOII, self.OIIIClIII, self.SIIIClIII, self.OIIIArVI, self.SIIIArVI, self.NIIClIII] :
          atr += lst

class Ion_Abun:

    def __init__(self, sizex=0, sizey=0):
        if sizex == 0 and sizey == 0:
            self.HeIa = []
            self.HeIb = []
            self.HeIIa = []
            self.HeIIb = []
            self.NI = []
            self.NIIa = []
            self.NIIb = []
            self.NIIc = []
            self.OIa = []
            self.OIb = []
            self.OIc = []
            self.OIIa = []
            self.OIIb = []
            self.OIIc = []
            self.OIId = []
            self.OIIIa = []
            self.OIIIb = []
            self.OIIIc = []
            self.SIIa = []
            self.SIIb = []
            self.SIIIa = []
            self.SIIIb = []
            self.NeIIIa = []
            self.NeIIIb = []
            self.ArIII = []
            self.ArIVa = []
            self.ArIVb = []
            self.ClIIIa = []
            self.ClIIIb = []
        else:
            self.HeIa = np.zeros((sizex, sizey))
            self.HeIb = np.zeros((sizex, sizey))
            self.HeIIa = np.zeros((sizex, sizey))
            self.HeIIb = np.zeros((sizex, sizey))
            self.NI = np.zeros((sizex, sizey))
            self.NIIa = np.zeros((sizex, sizey))
            self.NIIb = np.zeros((sizex, sizey))
            self.NIIc = np.zeros((sizex, sizey))
            self.OIa = np.zeros((sizex, sizey))
            self.OIb = np.zeros((sizex, sizey))
            self.OIc = np.zeros((sizex, sizey))
            self.OIIa = np.zeros((sizex, sizey))
            self.OIIb = np.zeros((sizex, sizey))
            self.OIIc = np.zeros((sizex, sizey))
            self.OIId = np.zeros((sizex, sizey))
            self.OIIIa = np.zeros((sizex, sizey))
            self.OIIIb = np.zeros((sizex, sizey))
            self.OIIIc = np.zeros((sizex, sizey))
            self.SIIa = np.zeros((sizex, sizey))
            self.SIIb = np.zeros((sizex, sizey))
            self.SIIIa = np.zeros((sizex, sizey))
            self.SIIIb = np.zeros((sizex, sizey))
            self.NeIIIa = np.zeros((sizex, sizey))
            self.NeIIIb = np.zeros((sizex, sizey))
            self.ArIII = np.zeros((sizex, sizey))
            self.ArIVa = np.zeros((sizex, sizey))
            self.ArIVb = np.zeros((sizex, sizey))
            self.ClIIIa = np.zeros((sizex, sizey))
            self.ClIIIb = np.zeros((sizex, sizey))


class Elem_Abun_KB:

    def __init__(self, sizex=0, sizey=0):
        if sizex == 0 and sizey == 0:
            self.He = []
            self.N = []
            self.O = []
            self.S = []
            self.Ne = []
            self.Ar = []
            self.Cl = []
        else:
            self.He = np.zeros((sizex, sizey))
            self.N = np.zeros((sizex, sizey))
            self.O = np.zeros((sizex, sizey))
            self.S = np.zeros((sizex, sizey))
            self.Ne = np.zeros((sizex, sizey))
            self.Ar = np.zeros((sizex, sizey))
            self.Cl = np.zeros((sizex, sizey))


class Elem_Abun_KB_Ratio:

    def __init__(self, sizex=0, sizey=0):
        if sizex == 0 and sizey == 0:
            self.NO = []
            self.SO = []
            self.NeO = []
            self.ArO = []
            self.ClO = []
        else:
            self.NO = np.zeros((sizex, sizey))
            self.SO = np.zeros((sizex, sizey))
            self.NeO = np.zeros((sizex, sizey))
            self.ArO = np.zeros((sizex, sizey))
            self.ClO = np.zeros((sizex, sizey))


class Ratio_Angles:

    def __init__(self):
        self.HeIa5876_Ha = []
        self.HeIb6678_Ha = []
        self.HeIIa4686_Hb = []
        self.HeIIb5412_Hb = []
        self.HeIa5876_HeIIa4686 = []
        self.HeIa5876_HeIIb5412 = []
        self.NIIc6583_Ha = []
        self.NIIbc654884_Ha = []
        self.NIIbc654884_NIIa5755 = []
        self.NIIbc654884_OIIIbc_495907 = []
        self.NI5199_Hb = []
        self.SIIab671631_Ha = []
        self.SIIa6716_SIIb6731 = []
        self.SIIab671631_SIIIab_631269 = []
        self.OIbc630063_Ha = []
        self.OIbc630063_OIa5577 = []
        self.OIbc630063_OIIIbc495907 = []
        self.OIbc630063_OIIcd732030 = []
        self.OIbc630063_OIIab372729 = []
        self.OIb6300_Ha = []
        self.OIIIc5007_Hb = []
        self.OIIIbc495907_Hb = []
        self.OIIIbc45907_OIIIa4363 = []
        self.OIIab372729_Hb = []
        self.OIIab372729_OIIIbc495907 = []
        self.OIIcd732030_OIIIbc495907 = []
        self.ArIVa4712_ArIVb4740 = []
        self.ArIVab471240_Hb = []
        self.NeIIIab386867_Hb = []
        self.ClIIIa5517_ClIIIb5538 = []
        self.ClIIIab551738_Hb = []
        self.CI8727_Ha = []
        self.CII6461_Ha = []

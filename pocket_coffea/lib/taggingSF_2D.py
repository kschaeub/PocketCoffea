import numpy as np
import awkward as ak
import correctionlib
import functools
import gzip
import cloudpickle
import correctionlib.schemav2 import Correction, CorrectionSet
import pocket_coffea.lib.correction_cache import load_correction



def get_scores(jets):
    """
    Computes new tagger scores based on UParTAK4 for the 2D pseudocontinuos SFs.
    """

    #get UParTAK4 scores
    upart_B = jets.btagUParTAK4B
    upart_SvUDG = jets.btagUParTAK4SvUDG
    upart_CvL = jets.btagUParTAK4CvL
    upart_CvB = jets.btagUParTAK4CvB

    #compute probabailites
    prob_udg = jets.UParTAK4UDG
    prob_s = ak.where(upart_B < 0, -1, upart_SvUDG * prob_udg / (1. - upart_SvUDG))
    prob_c = ak.where(upart_B < 0, -1, upart_CvL * (prob_s + prob_udg) / (1. - upart_CvL))
    prob_b = ak.where(upart_B < 0, -1, (1 - upart_CvB) * prob_c / (upart_CvB))


    #define new scores
    upart_BvsC = ak.where( upart_B < 0, -1, (1 - upart_CvB))
    upart_HFvsLF = ak.where( upart_B < 0, -1, (prob_b + prob_c)/ ( prob_b + prob_c + prob_s + prob_udg))


    jets["btagUParTAK4BvsC"] = upart_BvsC
    jets["btagUParTAK4HFvsLF"] = upart_HFvsLF

    return jets



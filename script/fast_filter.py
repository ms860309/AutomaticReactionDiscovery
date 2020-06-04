from makeit.utilities.fastfilter_utilities import Highway_self, pos_ct, true_pos, real_pos, set_keras_backend
from makeit.utilities.fingerprinting import create_rxn_Morgan2FP_separately
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from makeit.interfaces.scorer import Scorer
import numpy as np
import csv
from pymongo import MongoClient
from tqdm import tqdm
from keras.models import load_model
from keras import backend as K
import makeit.global_config as gc
from makeit.utilities.io.logger import MyLogger
import os
fast_filter_loc = 'fast_filter'


class FastFilterScorer(Scorer):
    def __init__(self):
        self.model = None

    def set_keras_backend(self, backend):
        if K.backend() != backend:
            os.environ['KERAS_BACKEND'] = backend
            reload(K)
            assert K.backend() == backend

    def load(self, model_path):
        MyLogger.print_and_log('Starting to load fast filter', fast_filter_loc)
        self.model = load_model(model_path, custom_objects={
                                'Highway_self': Highway_self, 'pos_ct': pos_ct, 'true_pos': true_pos, 'real_pos': real_pos})
        self.model._make_predict_function()
        MyLogger.print_and_log('Done loading fast filter', fast_filter_loc)

    def evaluate(self, reactant_smiles, target, **kwargs):
        # Strip chirality
        # rmol = Chem.MolFromSmiles(reactant_smiles)
        # pmol = Chem.MolFromSmiles(target)
        # reactant_smiles = Chem.MolToSmiles(rmol, False)
        # target = Chem.MolToSmiles(pmol, False)

        [pfp, rfp] = create_rxn_Morgan2FP_separately(
            reactant_smiles, target, rxnfpsize=2048, pfpsize=2048, useFeatures=False)
        pfp = np.asarray(pfp, dtype='float32')
        rfp = np.asarray(rfp, dtype='float32')
        rxnfp = pfp - rfp

        score = self.model.predict(
            [pfp.reshape(1, 2048), rxnfp.reshape(1, 2048)])
        outcome = {'smiles': target,
                   'template_ids': [],
                   'num_examples': 0
                   }
        all_outcomes = []
        all_outcomes.append([{'rank': 1.0,
                              'outcome': outcome,
                              'score': float(score[0][0]),
                              'prob': float(score[0][0]),
                              }])
        return all_outcomes

    def filter_with_threshold(self, reactant_smiles, target, threshold):
        [pfp, rfp] = create_rxn_Morgan2FP_separately(
            reactant_smiles, target, rxnfpsize=2048, pfpsize=2048, useFeatures=False)
        pfp = np.asarray(pfp, dtype='float32')
        rfp = np.asarray(rfp, dtype='float32')
        rxnfp = pfp - rfp

        score = self.model.predict([pfp.reshape(1, 2048), rxnfp.reshape(1, 2048)])
        filter_flag = (score > threshold)
        return filter_flag, float(score)


if __name__ == "__main__":

    ff = FastFilterScorer()
    ff.load(model_path=gc.FAST_FILTER_MODEL['trained_model_path'])
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'O=C[C@H](/C=C(/CO)\O)O.O')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'O[C@@H]1C(=O)[C@H]([C@@H]([C@H]1O)O)O.[H][H]')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'OC[C@@H](CC(C=O)(O)O)O')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'OCO[C@@H]([C@H](CO)O)C=O')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'OC[C@@H]1OC(=O)[C@H]([C@H]1O)O.[H][H]')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'O=C[C@H](C[C@H](C(O)O)O)O')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'OCC[C@@H](C(C=O)(O)O)O')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'O[C@H]1CC(=O)[C@H]([C@@H]1O)O.O')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'O[C@@H]([C@H](CO)O)OCC=O')
    print(score)
    score = ff.evaluate('OCC(C(C(C=O)O)O)O', 'O=CO[C@H]([C@@H](CO)O)CO')
    print(score)
    
    """
    flag, sco = ff.filter_with_threshold('CCO.CC(=O)O', 'CCOC(=O)C', 0.75)
    print(flag)
    print(sco)
    """
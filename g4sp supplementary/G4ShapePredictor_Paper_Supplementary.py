import pandas as pd
import os
import logomaker
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

from matplotlib.ticker import MultipleLocator
from tqdm import tqdm
from matplotlib.ticker import MultipleLocator
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import accuracy_score, roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from scipy.signal import find_peaks

#%%

def load_pickle(filepath): 
    return pickle.load(open(filepath, 'rb'))

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
df = pd.read_excel(os.path.join(script_dir, 'G4 Dataset.xlsx'), index_col=0)
threshold_params = load_pickle(os.path.join(script_dir, 'optimized_threshold.pkl'))
ova_threshold_params = load_pickle(os.path.join(script_dir, 'ova_optimized_threshold.pkl'))

CBC_param_grid= {
    'iterations': [363],
    'boosting_type': ['Plain'],
    'depth': [6],
}

EXT_param_grid = {
    'n_estimators': [900], 
    'criterion': ['entropy'],
    'max_depth': [14],
    'min_samples_split': [5],
    'min_samples_leaf': [1],
    'max_features': ['sqrt'],
    'max_leaf_nodes': [None],
    'bootstrap': [True],
    'class_weight': ['balanced_subsample'],
    'max_samples': [None]
}

RF_param_grid = {
    'n_estimators': [182],
    'criterion': ['gini'],
    'max_depth': [15],
    'min_samples_split': [2],
    'min_samples_leaf': [1],
    'min_weight_fraction_leaf': [0.0],
    'max_features': ['sqrt'],
    'max_leaf_nodes': [None],
    'bootstrap': [True],
    'class_weight': ['balanced'],
}

XGB_param_grid = {
    'booster': ['gbtree'],
    'eta': [0.1],
    'gamma': [0],
    'max_depth': [6],
    'min_child_weight': [1],
    'subsample': [0.7],
    'tree_method': ['approx'],
    'process_type': ['default'],
}

GBC_param_grid = {
    'learning_rate': [0.239],
    'n_estimators': [210],  
    'criterion': ['friedman_mse'],
    'min_samples_split': [12],  
    'min_samples_leaf': [1],
    'max_depth': [10], 
    'min_impurity_decrease': [0], 
    'max_features': ['sqrt'],
    'max_leaf_nodes': [70], 
    'validation_fraction': [0.1],  
    'n_iter_no_change': [None],  
    'tol': [1e-4],   
}

LGBM_param_grid = {
    'n_estimators': [342],
    'boosting_type': ['gbdt'],
    'learning_rate': [0.069],
    'num_leaves': [28],
    'tree_learner': ['serial'],
    'application':['multiclassova'],
    'num_class':[1]
}

param_grids = [EXT_param_grid, RF_param_grid, XGB_param_grid, GBC_param_grid, LGBM_param_grid, CBC_param_grid]

for param_grid in param_grids:
    for key, value in param_grid.items():
        if isinstance(value, list):
            param_grid[key] = value[0]
            
def flatten(c):
    return [a for b in c for a in b]

def pad(sequence, pad_length=100):
    while len(sequence) < pad_length:
        sequence = 'N' + sequence + 'N'
    while len(sequence) > pad_length:
        sequence = sequence[:-1]
    return sequence

def close(value, dct):
    """
    Find the key in the dictionary that is closest to the provided value.
    
    Parameters:
    - value: The value to compare against.
    - dct: The dictionary whose keys are to be compared.
    
    Returns:
    - The key from the dictionary that is closest to the provided value.
    """
    return min(dct.keys(), key=lambda x: abs(x - value))

def seqs_converter(input_seqs, pad_length=100): # converts LIST of ATCGU = 1,2,3,4,5 and shape (,100)
    if type(input_seqs) == str:
        input_seqs = [input_seqs]
    seqs = []
    seqsdict = {'A':1,'T':2,'C':3,'G':4,'U':5,'N':0,'Y':0}
    for i in input_seqs:
        temp = [seqsdict[x] for x in i]
        while len(temp) < pad_length:
            temp.insert(0,0)
            temp.append(0)
        while len(temp) > pad_length:
            del temp[-1]
        seqs.append(temp)
    return [np.array(x) for x in seqs]

def rseqs_converter(input_seqs): # converts LIST of 0,1,2,3,4 = N, A, T, C, G and shape (,100) back to ATCGN
    seqs_inverse_dict = {1:'A', 2:'T', 3:'C', 4:'G', 5:'U', 0:''}
    return [''.join(seqs_inverse_dict[x] for x in i) for i in input_seqs]

class G4Data:
    
    def __init__(self, df):
        self.df = df
        
        X, y = np.array(seqs_converter(list(df['Sequence']), pad_length=100)), np.asarray(df['Topology'])
        self.X = X
        self.y = y
        
        kfold = StratifiedKFold(random_state=42, n_splits=10, shuffle=True)
        self.kfold = kfold
        
        repeat_kfold = RepeatedStratifiedKFold(n_splits=10, n_repeats=100, random_state=42)
        self.repeat_kfold = repeat_kfold
        
        self.models = ['CatBoostClassifier (CBC)',
                       'ExtraTreesClassifier (EXT)', 
                       'GradientBoostingClassifier (GBC)',
                       'Light Gradient-boosted Machine (LGBM)',
                       'RandomForest (RF)',
                       'XGBoost (XGB)']
    
    def display_models(self):
        print('Models are:\n' + ',\n'.join(self.models))
        
    def PDB_data(self):
        ''' 
        Returns G4 data extracted from PDB: https://www.rcsb.org
        '''
        return self.df[self.df['Comments'].apply(lambda x: len(x) == 4)]
    
    def lit_data(self):
        ''' 
        Returns G4 data extracted from literature
        and NOT FOUND in PDB
        '''
        return self.df[self.df['Comments'].apply(lambda x: 'In-house CD' not in x and len(x) > 4)]
    
    def inhouse_data(self):
        '''
        Returns G4 data obtained from in-house CD
        experiments
        '''
        return self.df[self.df['Comments'].apply(lambda x: all(y in x for y in ['PQS', 'In-house CD']))]
    
    def parallel_G4(self):
        '''
        Returns all parallel G4 data from all sources
        '''
        return self.df[self.df['Topology'] == 0]
    
    def antiparallel_G4(self):
        '''
        Returns all antiparallel G4 data from all sources
        '''
        return self.df[self.df['Topology'] == 1]
    
    def hybrid_G4(self):
        '''
        Returns all hybrid G4 data from all sources
        '''
        return self.df[self.df['Topology'] == 2]
    
    def sequence_lengths(self):
        '''
        Returns a list of lengths of input G4 sequence
        '''
        df = self.df.copy()
        df['Length'] = df['Sequence'].apply(len)
        return df[['Length', 'Sequence', 'Topology']]
    
    def get_mean(self):
        '''
        Returns mean length of input G4 sequences
        '''
        df = self.df.copy()
        return df['Sequence'].apply(len).mean()
    
    def get_median(self, df):
        '''
        Returns median length of input G4 sequences
        '''
        return df['Sequence'].apply(len).median()
    
    def get_mode(self, df):
        '''
        Returns mode length of input G4 sequences
        '''
        return df['Sequence'].apply(len).mode()
    
    def get_std(self, df):
        '''
        Returns standard deviation of lengths of input G4 sequences
        '''
        return df['Sequence'].apply(len).std()
    
    def get_iqr(self, df):
        '''
        Returns interquartile range of lengths of input G4 sequences
        '''
        lengths = df['Sequence'].apply(len)
        Q1 = lengths.quantile(0.25)
        Q3 = lengths.quantile(0.75)
        return Q3 - Q1
    
    def frequency_plot(self, topology, show=True):
        topology_garden = {0: self.parallel_G4(),
                            1: self.antiparallel_G4(),
                            2: self.hybrid_G4()}
        titles = {0:'Parallel (4+0)',
                  1:'Antiparallel (2+2)',
                  2:'Hybrid (3+1)'}
        title = titles[topology]
        
        sequence_list = topology_garden[topology]['Sequence'].tolist()
        
        data = {'A': [0 for _ in range(100)],
                'T': [0 for _ in range(100)],
                'C': [0 for _ in range(100)],
                'G': [0 for _ in range(100)],
        }
        
        for sequence in sequence_list:
            sequence = pad(sequence)
            for ndx, nucleotide in enumerate(flatten(sequence)):
                if nucleotide !='N':
                    data[nucleotide][ndx] += 1
        
        df = pd.DataFrame(data).transpose()
        df_normalized = df.div(df.sum(axis=0), axis=1).fillna(0)
        df_plot = df_normalized.transpose().copy()
        df_plot.index += 1
        
        logomaker.Logo(df_plot, color_scheme='classic',
                                font_name='Arial')
        
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.xlim([35,65])
        plt.title(f'Logo Plot Representation of {title}')
        
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_tick_params(which='minor', size=5, direction='out', top=False)
        if show:
            plt.show()
        else:
            plt.clf()
        return df
    
    def frequency_plot_parallel(self):
        return self.frequency_plot(0)
    
    def frequency_plot_antiparallel(self):
        return self.frequency_plot(1)
        
    def frequency_plot_hybrid(self):
        return self.frequency_plot(2)
    
    def frequency_plot_linegraph(self, topology):
        df = self.frequency_plot(topology, show=False).transpose()
        color_garden = {'A': '#008000',
                        'T': '#ff0000',
                        'C': '#2727ff',
                        'G': '#ffa600'}
        df = df/max(np.array(df).flatten())
        
        plt.figure(dpi=300)
        for col in df.columns:
            plt.plot(range(1, 101), df[col], label=col, color=color_garden[col])
        plt.xlim(35, 65)
        plt.legend()
        plt.xlabel('Nucleotide Position', fontsize=12)
        plt.ylabel('Nucleotide Frequency', fontsize=12)

        ax = plt.gca()
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        
    def frequency_plot_linegraph_parallel(self):
       return self.frequency_plot_linegraph(0)
        
    def frequency_plot_linegraph_antiparallel(self):
        return self.frequency_plot_linegraph(1)
    
    def frequency_plot_linegraph_hybrid(self):
        return self.frequency_plot_linegraph(2)
        
    def model_garden(self):
        models = {
            'CBC': CatBoostClassifier(logging_level='Silent', random_state=42, **CBC_param_grid),
            'EXT': ExtraTreesClassifier(random_state=42, **EXT_param_grid, n_jobs=-1),
            'GBC': GradientBoostingClassifier(random_state=42, **GBC_param_grid), 
            'LGBM': LGBMClassifier(random_state=42, **LGBM_param_grid, n_jobs=-1), 
            'RF': RandomForestClassifier(random_state=42, **RF_param_grid, n_jobs=-1), 
            'XGB': XGBClassifier(random_state=42, **XGB_param_grid, n_jobs=-1)
            }
        return models
    
    def topology_garden(self):
        topologies = {0: 'parallel',
                      1: 'antiparallel',
                      2: 'hybrid'}
        return topologies
        
    def train(self, model_name):
        scores = []
        X, y = self.X.copy(), self.y.copy()
        model = self.model_garden()[model_name]
        print(f'Training {model_name} with K-Fold Cross-Validation of 10 folds, with each fold being a random shuffle of the dataset')
        
        for train_index, test_index in tqdm(self.kfold.split(X, y), desc='Training with 10 folds ...', total=10):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            scores.append(accuracy_score(y_test, y_pred))
        
        return {'mean accuracy': np.mean(scores),
                'std': np.std(scores)}
    
    def train_OvA(self, model_name):
        X, y = self.X.copy(), self.y.copy()
        classes = [0, 1, 2]
        
        scores = [[] for _ in classes]
        model = OneVsRestClassifier(self.model_garden()[model_name])
        print(f'Training {model_name} (OneVsRestClassifier) with K-Fold Cross-Validation of 10 folds, with each fold being a random shuffle of the dataset')
        for train_index, test_index in tqdm(self.kfold.split(X, y), desc='Training with 10 folds ...', total=10):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            
            model.fit(X_train, y_train)
            
            y_pred = model.predict(X_test)
            for class_index in range(len(classes)):
                y_true_class = (y_test == class_index).astype(int)
                y_pred_class = (y_pred == class_index).astype(int)
                accuracy = accuracy_score(y_true_class, y_pred_class)
                scores[class_index].append(accuracy)
        
        return {'model': model_name,
                'mean accuracy of parallel': np.mean(scores[0]),
                'std accuracy of parallel': np.std(scores[0]),
                'mean accuracy of antiparallel': np.mean(scores[1]),
                'std accuracy of antiparallel': np.std(scores[1]),
                'mean accuracy of hybrid': np.mean(scores[2]),
                'std accuracy of hybrid': np.std(scores[2])}
                
    def train_CBC(self):
        return self.train('CBC')
        
    def train_EXT(self):
        return self.train('EXT')
    
    def train_GBC(self):
        return self.train('GBC')
    
    def train_LGBM(self):
        return self.train('LGBM')
    
    def train_RF(self):
        return self.train('RF')
    
    def train_XGB(self):
        return self.train('XGB')
    
    def train_ova_CBC(self):
        return self.train_OvA('CBC')
    
    def train_ova_EXT(self):
        return self.train_OvA('EXT')
    
    def train_ova_GBC(self):
        return self.train_OvA('GBC')
    
    def train_ova_LGBM(self):
        return self.train_OvA('LGBM')
    
    def train_ova_RF(self):
        return self.train_OvA('RF')
    
    def train_ova_XGB(self):
        return self.train_OvA('XGB')
    
    def plot_auroc(self, model_name):
        
        X, y = self.X.copy(), self.y.copy()
        model_garden = {'CBC': CatBoostClassifier(logging_level='Silent', random_state=42, **CBC_param_grid),
                      'EXT': ExtraTreesClassifier(random_state=42, **EXT_param_grid, n_jobs=-1),
                      'GBC': GradientBoostingClassifier(random_state=42, **GBC_param_grid), 
                      'LGBM': LGBMClassifier(random_state=42, **LGBM_param_grid, n_jobs=-1), 
                      'RF': RandomForestClassifier(random_state=42, **RF_param_grid, n_jobs=-1), 
                      'XGB': XGBClassifier(random_state=42, **XGB_param_grid, n_jobs=-1)}
        
        model = model_garden[model_name]

        classes = [0, 1, 2]
        y_bin = label_binarize(y, classes=classes)
        n_classes = y_bin.shape[1]
        
        fpr, tpr, roc_auc = dict(), dict(), dict()

        classifier = OneVsRestClassifier(model)
        
        for train_index, val_index in tqdm(self.repeat_kfold.split(X, y)):
            X_train, X_val = X[train_index], X[val_index]
            y_train, y_val = y_bin[train_index], y_bin[val_index]
        
            classifier.fit(X_train, y_train)
            y_score = classifier.predict_proba(X_val)
        
            for i in range(n_classes):
                fpr[i], tpr[i], _ = roc_curve(y_val[:, i], y_score[:, i])
                roc_auc[i] = auc(fpr[i], tpr[i])
        
        colors = ['#292FFF', '#34FF34', '#FF1A1A']
        plt.figure(dpi = 300)
        lines = []
        labels = []
        classes = ['(4+0)', '(2+2)', '(3+1)']
        for i, color in zip(range(n_classes), colors):
            l, = plt.plot(fpr[i], tpr[i], color=color, lw=2)
            lines.append(l)
            labels.append('{0} (area = {1:0.2f})'.format(classes[i], roc_auc[i]))
            
        plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
        
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'Receiver Operating Characteristic ({model_name})')
        plt.legend(lines, labels, loc='lower right', prop=dict(size=10))
        plt.tight_layout()
        plt.savefig(rf'C:\Desktop\Kitchen Sink\{model_name}_AUROC.pdf', dpi=300, bbox_inches='tight')
        plt.show()
        
    def AUROC_CBC(self):
        self.plot_auroc('CBC')
    
    def AUROC_EXT(self):
        self.plot_auroc('EXT')
    
    def AUROC_GBC(self):
        self.plot_auroc('GBC')
    
    def AUROC_LGBM(self):
        self.plot_auroc('LGBM')
    
    def AUROC_RF(self):
        self.plot_auroc('RF')
    
    def AUROC_XGB(self):
        self.plot_auroc('XGB')
    
    def plot_PR(self, model_name):
        
        X, y = self.X.copy(), self.y.copy()
        model_garden = {'CBC': CatBoostClassifier(logging_level='Silent', random_state=42, **CBC_param_grid),
                      'EXT': ExtraTreesClassifier(random_state=42, **EXT_param_grid, n_jobs=-1),
                      'GBC': GradientBoostingClassifier(random_state=42, **GBC_param_grid), 
                      'LGBM': LGBMClassifier(random_state=42, **LGBM_param_grid, n_jobs=-1), 
                      'RF': RandomForestClassifier(random_state=42, **RF_param_grid, n_jobs=-1), 
                      'XGB': XGBClassifier(random_state=42, **XGB_param_grid, n_jobs=-1)}
        
        model = model_garden[model_name]
        
        classes = [0, 1, 2]
        
        y_bin = label_binarize(y, classes=classes)
        n_classes = y_bin.shape[1]
        
        precision, recall, threshold, average_precision = dict(), dict(), dict(), dict()
        
        classifier = OneVsRestClassifier(model)
        
        for train_index, val_index in tqdm(self.repeat_kfold.split(X, y)):
            X_train, X_val = X[train_index], X[val_index]
            y_train, y_val = y_bin[train_index], y_bin[val_index]
        
            classifier.fit(X_train, y_train)
            y_score = classifier.predict_proba(X_val)
        
            for i in range(n_classes):
                precision[i], recall[i], threshold[i] = precision_recall_curve(y_val[:, i], y_score[:, i])
                precision[i], recall[i] = precision[i][:-1], recall[i][:-1] # drop the last element of PR
                average_precision[i] = average_precision_score(y_val[:, i], y_score[:, i])
        
        colors = ['#292FFF', '#34FF34', '#FF1A1A']
        plt.figure(dpi=300)
        lines = []
        labels = []
        classes = ['(4+0)', '(2+2)', '(3+1)']
        for i, color in zip(range(n_classes), colors):
            l, = plt.plot(recall[i], precision[i], color=color, lw=2)
            lines.append(l)
            labels.append('{0} (area = {1:0.2f})'.format(classes[i], average_precision[i]))
        
        for i in range(n_classes):
            no_skill = sum(y == i) / len(y)
            plt.plot([0, 1], [no_skill, no_skill], linestyle='--', color=colors[i])
        
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(f'Precision-Recall curve ({model_name})')
        plt.tight_layout()
        plt.legend(lines, labels, loc='best', prop=dict(size=10))
        plt.tight_layout()
        plt.savefig(rf'C:\Desktop\Kitchen Sink\{model_name}_AUPR.pdf', dpi=300, bbox_inches='tight')
        plt.show()
        
    def PR_CBC(self):
        self.plot_PR('CBC')
    
    def PR_EXT(self):
        self.plot_PR('EXT')
    
    def PR_GBC(self):
        self.plot_PR('GBC')
    
    def PR_LGBM(self):
        self.plot_PR('LGBM')
    
    def PR_RF(self):
        self.plot_PR('RF')
    
    def PR_XGB(self):
        self.plot_PR('XGB')
        
    def train_optimized_threshold(self, model_name, threshold=0.9):
        scores, thresholds = [], []
        X, y = self.X.copy(), self.y.copy()
        model = self.model_garden()[model_name]
        for topology in [0, 1, 2]: # 0: parallel, 1: antiparallel, 2: hybrid
            desired_precision = close(threshold, threshold_params[model_name][topology])
            thresholds.append(threshold_params[model_name][topology][desired_precision]['threshold'])
        
        print(f'Optimized thresholds of {thresholds[0]},{thresholds[1]},{thresholds[2]} for parallel,antiparallel,hybrid respectively have been selected.')
        for train_index, test_index in tqdm(self.repeat_kfold.split(X, y), desc='Training with 10 folds ...', total=10):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            
            model.fit(X_train, y_train)
            y_pred_proba = model.predict_proba(X_test)
            
            filtered_list = {'true':[], 'pred':[]}
            for idx, x in enumerate(y_pred_proba):
                if np.sum([x[i] > thresholds[i] for i in range(len(thresholds))]) == 1:
                    filtered_list['true'].append(y_test[idx])
                    filtered_list['pred'].append(np.argmax(y_pred_proba[idx]))
                    
            scores.append(accuracy_score(filtered_list['true'], filtered_list['pred']))
        return {'mean accuracy': np.mean(scores),
                'std': np.std(scores),
                'model': model_name}
    
    def train_optimized_threshold_OvA(self, model_name, threshold=0.9):
        scores, thresholds = {0:[], 1:[], 2:[]}, []
        X, y = self.X.copy(), self.y.copy()
        base_model = self.model_garden()[model_name]
        ova_model = OneVsRestClassifier(base_model)
        for topology in [0, 1, 2]: # 0: parallel, 1: antiparallel, 2: hybrid
            desired_precision = close(threshold, ova_threshold_params[model_name][topology])
            thresholds.append(ova_threshold_params[model_name][topology][desired_precision]['threshold'])
        
        print(f'Optimized thresholds of {thresholds[0]},{thresholds[1]},{thresholds[2]} for parallel, antiparallel, hybrid respectively have been selected.')
        for train_index, test_index in tqdm(self.repeat_kfold.split(X, y), desc='Training with 10 folds ...', total=10):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            
            ova_model.fit(X_train, y_train)
            
            temp_scores = {0:[], 1:[], 2:[]}
            for topology in [0, 1, 2]:
                binary_classifier = ova_model.estimators_[topology]
                positive_class_probability = [x[1] for x in binary_classifier.predict_proba(X_test)]
                
                for idx, x in enumerate(positive_class_probability):
                    if positive_class_probability[idx] > thresholds[topology]:
                        if topology == y_test[idx]:
                            temp_scores[topology].append(True)
                        else:
                            temp_scores[topology].append(False)
            
            for topology in [0, 1, 2]:
                if temp_scores[topology]:
                    scores[topology].append(np.mean(temp_scores[topology]))
        
        mean_accuracy = {topology: np.mean(scores[topology]) for topology in scores}
        std_accuracy = {topology: np.std(scores[topology]) for topology in scores}
        return {'mean accuracy': mean_accuracy, 
                'std accuracy': std_accuracy,
                'model': model_name}
    
    def tt(self, model_name, threshold=0.9):
        thresholds = []
        base_model = self.model_garden()[model_name]
        ova_model = OneVsRestClassifier(base_model)
    
        for topology in [0, 1, 2]:  # 0: parallel, 1: antiparallel, 2: hybrid
            desired_precision = close(threshold, threshold_params[model_name][topology])
            thresholds.append(threshold_params[model_name][topology][desired_precision]['threshold'])
    
        print(f'Optimized thresholds of {thresholds[0]},{thresholds[1]},{thresholds[2]} for parallel, antiparallel, hybrid respectively have been selected.')
    
        scores = {0: [], 1: [], 2: []}  # Dictionary to hold scores for each class
    
        for i in tqdm(range(len(self.X)), desc=f'Training {model_name} using optimized thresholds ...'):
            temp_test_X, temp_test_y = self.X[i, :].reshape(1, -1), self.y[i]
            temp_train_X, temp_train_y = np.delete(self.X, i, axis=0), np.delete(self.y, i)
            ova_model.fit(temp_train_X, temp_train_y)
    
            # Extract binary classifiers
            for topology in [0, 1, 2]:
                binary_classifier = ova_model.estimators_[topology]
                prediction_probability = binary_classifier.predict_proba(temp_test_X)[0][1]  # Probability of class 1
                
                if prediction_probability > thresholds[topology]:
                    scores[topology].append(temp_test_y == topology)
                    
                else:
                    scores[topology].append(not temp_test_y != topology)
                
        mean_accuracy = {topology: np.mean(scores[topology]) for topology in scores}
    
        return {'mean accuracy': mean_accuracy, 'model': model_name}
        
    def train_optimized_CBC(self, threshold=0.9):
        return self.train_optimized_threshold('CBC', threshold=threshold)
    
    def train_optimized_EXT(self, threshold=0.9):
        return self.train_optimized_threshold('EXT', threshold=threshold)
    
    def train_optimized_GBC(self, threshold=0.9):
        return self.train_optimized_threshold('GBC', threshold=threshold)
    
    def train_optimized_LGBM(self, threshold=0.9):
        return self.train_optimized_threshold('LGBM', threshold=threshold)
    
    def train_optimized_RF(self, threshold=0.9):
        return self.train_optimized_threshold('RF', threshold=threshold)
    
    def train_optimized_XGB(self, threshold=0.9):
        return self.train_optimized_threshold('XGB', threshold=threshold)
    
    def get_metrics_vs_threshold(self, model_name):
        
        X, y = self.X.copy(), self.y.copy()
        model = self.model_garden()[model_name]
        classes = [0, 1, 2]
        metrics = []
        for topology in classes:
            y_vals, y_scores = [], []
            for train_index, val_index in tqdm(self.kfold.split(X, y), desc=f'Getting metrics for {self.topology_garden()[topology]} ...', total=self.kfold.n_splits):
                X_train, X_val = X[train_index], X[val_index]
                y_train, y_val = y[train_index], y[val_index]
                
                model.fit(X_train, y_train)
                predictions = model.predict_proba(X_val)
                y_vals.append(y_val)
                y_scores.append(predictions)
                
            precision, recall, thresholds = precision_recall_curve(np.concatenate(y_vals) == topology, np.concatenate(y_scores,axis=0)[:, topology])
            df = pd.DataFrame([thresholds, precision, recall]).transpose()
            df.columns = [f'Threshold_{topology}', 
                          f'Precision_{topology}', 
                          f'Recall_{topology}']
            metrics.append(df)
        return pd.concat(metrics, axis=1)
    
    def get_metrics_vs_threshold_OvA(self, model_name):
        
        X, y = self.X.copy(), self.y.copy()
        model = self.model_garden()[model_name]
        ova_model = OneVsRestClassifier(model)
        classes = [0, 1, 2]
        metrics = []
        for topology in classes:
            y_tests, y_scores = [], []
            for train_index, test_index in tqdm(self.kfold.split(X, y), desc=f'Getting metrics for {self.topology_garden()[topology]} ...', total=self.kfold.n_splits):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                
                ova_model.fit(X_train, y_train)
                predictions = ova_model.estimators_[topology].predict_proba(X_test)
                y_tests.append(y_test)
                y_scores.append(predictions)
                
            precision, recall, thresholds = precision_recall_curve(np.concatenate(y_tests) == topology, np.concatenate(y_scores)[:, 1])
            df = pd.DataFrame([thresholds, precision, recall]).transpose()
            df.columns = [f'Threshold_{topology}', 
                          f'Precision_{topology}', 
                          f'Recall_{topology}']
            metrics.append(df)
        return pd.concat(metrics, axis=1)
    
    def get_metrics_vs_threshold_CBC(self):
        return self.get_metrics_vs_threshold('CBC')
    
    def get_metrics_vs_threshold_EXT(self):
        return self.get_metrics_vs_threshold('EXT')
        
    def get_metrics_vs_threshold_GBC(self):
        return self.get_metrics_vs_threshold('GBC')
    
    def get_metrics_vs_threshold_LGBM(self):
        return self.get_metrics_vs_threshold('LGBM')
    
    def get_metrics_vs_threshold_RF(self):
        return self.get_metrics_vs_threshold('RF')
    
    def get_metrics_vs_threshold_XGB(self):
        return self.get_metrics_vs_threshold('XGB')
        
    def get_feature_importances(self, model_name):
        
        X = self.X.copy()
        model = self.model_garden()[model_name]
        classes = [0, 1, 2]
        fold_importances = [[] for _ in range(len(classes))]
        rkfold = RepeatedStratifiedKFold(n_splits=10, n_repeats=10, random_state=42)
        for topology in classes:
            y = self.y.copy()
            y[y==topology] = -1 # placeholder
            y[y!=-1] = 0
            y[y==-1] = 1
            for train_idx, test_idx in tqdm(rkfold.split(X, y), total=100, desc=f'Obtaining feature importances for {model_name} ...'):
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]
                
                model.fit(X_train, y_train)
                perm_impt = permutation_importance(model, X_test, y_test, scoring='accuracy', n_repeats=10, random_state=42)
                fold_importances[topology].append(perm_impt.importances_mean)
        return fold_importances
            
    def identify_peaks(self, feature_importance, model_name, r=20):
        
        def plot_peaks(x, y, top=0, title= 'Detecting Peaks using scipy.signal.find_peaks',
               distance=20):
                """
                Plot data and highlight peaks using scipy's find_peaks.
            
                Parameters:
                - x: array-like, data for the x-axis
                - y: array-like, data for the y-axis
                - distance: int, minimum number of samples separating peaks
                """
                topology_color = {0:'#0000FF', 1:'#00FF00', 2:'#FF0000'}
                
                x, y = np.array(x), np.array(y)
                # Detect peaks
                peaks, _ = find_peaks(y, distance=distance)
            
                # Plot data and peaks
                plt.figure(figsize=(10, 6))
                plt.plot(x, y, label=self.topology_garden()[top], color = topology_color[top])
                plt.plot(x[peaks], y[peaks], "o", label='Peaks', color='#000000')
                plt.title(title)

                for peak in peaks:
                    plt.annotate(str(x[peak]), (x[peak], y[peak]), textcoords="offset points", xytext=(0,5), ha='center')
                
                plt.xlim(30, 70)
                plt.ylim(-0.005, 0.03)
                plt.legend()
                plt.show()
                
                
                return peaks
        
        classes = [0, 1, 2]
        feature_importance = [np.mean(x, axis=0) for x in feature_importance]
        peaks = []
        
        for topology in classes:
            peaks.append(plot_peaks(range(1, len(feature_importance[topology])+1), feature_importance[topology], top=topology, title=f'{model_name}', distance=r).tolist())
        
        return {self.topology_garden[c]:peaks[c] for c in classes}
    
    def feature_importances_CBC(self, distance=10):
        return self.identify_peaks(self.get_feature_importances('CBC'), 'CBC', r=distance)
    
    def feature_importances_EXT(self, distance=10):
        return self.identify_peaks(self.get_feature_importances('EXT'), 'EXT', r=distance)
    
    def feature_importances_GBC(self, distance=10):
        return self.identify_peaks(self.get_feature_importances('GBC'), 'GBC', r=distance)
    
    def feature_importances_LGBM(self, distance=10):
        return self.identify_peaks(self.get_feature_importances('LGBM'), 'LGBM', r=distance)
    
    def feature_importances_RF(self, distance=10):
        return self.identify_peaks(self.get_feature_importances('RF'), 'RF', r=distance)
    
    def feature_importances_XGB(self, distance=10):
        return self.identify_peaks(self.get_feature_importances('XGB'), 'XGB', r=distance)
    
    def permutate_sequences(self, length, peak): # only one peak at a time

        def centered_indices(p=peak, l=length):
            half_length = l // 2
            if l % 2 == 0:  # even length
                start = p - half_length
                end = p + half_length
            else:           # odd length
                start = p - half_length
                end = p + half_length + 1
            return list(range(start, end))
    

        X = self.X.copy()
        permutations = list(np.array(x) for x in tqdm(itertools.product([1,2,3,4], repeat=length), total = 4**length, desc='Permutating Sequences'))
        indices = centered_indices(peak, length)
    
        # Get a boolean array where each row is True if any value in that row of the subarray is zero
        rows_with_zeros = np.any(X[:, indices] == 0, axis=1)
        
        X = X[~rows_with_zeros]
        
        output_X = []
        for p in permutations:
            X_temp = X.copy()
            for i in range(X_temp.shape[0]):
                X_temp[i][indices] = p
            output_X.append(X_temp.copy())
        return (output_X, rseqs_converter(permutations))

    def heatmap_preparation_feature_importances(self, topology, peak, length):
        permutation_group = {}
        X, y = self.X.copy(), self.y.copy()
        X_permutations, permutations = self.permutate_sequences(length=length, peak=peak)
        for model_name, model in self.model_garden().items():
            
            trained_model = model.fit(X, y)
            
            for idx, x in enumerate(tqdm(X_permutations, desc = model_name)):
                prediction = trained_model.predict_proba(x)
                result = np.round(np.median(prediction, axis=0), 4).tolist()[topology]
                if permutations[idx] not in permutation_group:
                    permutation_group[permutations[idx]] = {}
                permutation_group[permutations[idx]][f'median_{model_name}'] = result
            
        return pd.DataFrame(permutation_group).transpose()
    
    def plot_heatmap(self, df, topology):
        length = len(df.index[0])
        df = df.mean(axis=1)
        heatmap_dict, indexes = {}, []
        topology_color = {0:'Blues', 1:'Greens', 2:'Reds'}
        
        for key, value in dict(df).items():
            prefix, suffix = key[:int(np.ceil(length/2))], key[-1*int(np.floor(length/2)):]

            indexes.append(suffix) if suffix not in indexes else None
            
            if prefix not in heatmap_dict:
                heatmap_dict[prefix] = []
            heatmap_dict[prefix].append(value)
        df = pd.DataFrame(heatmap_dict, index=indexes)
        df_normalized = (df - df.min().min()) / (df.max().max() - df.min().min())

        fig, ax = plt.subplots(figsize=(8, 12), dpi=300)
        
        ax = plt.gca()
        
        sns.heatmap(df_normalized.T, cbar=True, annot=False, cmap = (sns.color_palette(topology_color[topology], as_cmap=True)), yticklabels=True)

        ax.tick_params(axis='both', which='both', length=0)
        labels = [item.get_text()[:-1] + r'$\mathbf{' + item.get_text()[-1] + '}$' for item in ax.get_yticklabels()]
        ax.set_yticklabels(labels, fontsize=12, rotation=0)
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, rotation=0)
        
        plt.show()
        
        return df.transpose()     
        
    def get_heatmap(self, topology, peak, length):
        data = self.plot_heatmap(self.heatmap_preparation_feature_importances(topology=topology, 
                                                                              peak=peak,
                                                                              length=length), 
                                 topology=topology)
        return data
    

g = G4Data(df)



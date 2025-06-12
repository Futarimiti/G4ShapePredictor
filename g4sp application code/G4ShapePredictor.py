import PySimpleGUI as sg
from Bio import SeqIO
import pickle
import pandas as pd
import numpy as np
import os
import sys

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
# FIXME: catboost unspported on darwin
# from catboost import CatBoostClassifier
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier
#%%

# Get the directory of the .exe file
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
else:
    application_path = os.path.dirname(os.path.abspath(__file__))
    

def validate_sequence(sequence):
    return all(character in 'ATCGN' for character in sequence)

def seqs_converter(sequence_list, pad_length=100):
    
    seqsdict = {'A': 1, 'T': 2, 'C': 3, 'G': 4, 'N': 0}
    
    if isinstance(sequence_list, str):
        sequence_list = [sequence_list]
    
    def convert_sequence(sequence):
        if not validate_sequence(sequence):
            return 'Invalid'
        
        sequence_letters = [seqsdict[character] for character in sequence]
        
        # Add padding to the beginning and end of the sequence until it reaches the desired length
        while len(sequence_letters) < pad_length:
            sequence_letters.insert(0, 0)
            sequence_letters.append(0)
        
        # If the sequence is longer than the desired length, trim it from the end
        while len(sequence_letters) > pad_length:
            sequence_letters.pop()
        
        return np.array(sequence_letters)

    return [convert_sequence(sequence) for sequence in sequence_list]

def G4Hunter(line):
    scores = []
    length = len(line)
    i = 0
    
    while i < length:
        char = line[i]
        if char in ['G', 'g']:
            run_length = 1
            while i + run_length < length and line[i + run_length] in ['G', 'g']:
                run_length += 1
            scores.extend([run_length if run_length <= 4 else 4] * run_length)
            i += run_length

        elif char in ['C', 'c']:
            run_length = 1
            while i + run_length < length and line[i + run_length] in ['C', 'c']:
                run_length += 1
            scores.extend([-run_length if run_length <= 4 else -4] * run_length)
            i += run_length

        else:
            scores.append(0)
            i += 1

    return sum(scores) / length

def top_dict(inverted=False, ovr=None):
    topdict = {0: 'Parallel (4+0)', 1: 'Antiparallel (2+2)', 2: 'Hybrid (3+1)', -1: 'Mixed', -2: 'Invalid'}

    if isinstance(ovr, int):
        return {0: topdict[ovr], 1: f'Not {topdict[ovr]}', -2: 'Invalid'}

    if inverted:
        return topdict

    return {v: k for k, v in topdict.items()}
    
def is_valid_fasta(filename):
    try:
        valid_chars = {'A', 'T', 'C', 'G', 'N'}
        return all(set(seq.seq.upper()) <= valid_chars for seq in SeqIO.parse(filename, 'fasta'))
    except Exception as e:
        print(f'Error when reading file: {e}')
        return False
    
def load_model(filepath): # opens ML model using pickle
    return pickle.load(open(filepath, 'rb'))

def predict_G4(list_of_sequences, model, threshold='', ovr = False):
    
    def threshold_test(list_of_proba, threshold):
        return any([x>y for x,y in zip(list_of_proba, threshold)])
    
    # Handle 'Invalid characters'
    valid_indices = [i for i, seq in enumerate(list_of_sequences) if str(seq) != 'Invalid']
    valid_sequences = [list_of_sequences[i] for i in valid_indices]

    results = [-2] * len(list_of_sequences)  # Initialize results with all 'Invalid'

    if valid_sequences:  # Check if there are any valid sequences
        if not ovr:
            if bool(threshold):
                probabilities = model.predict_proba(valid_sequences)
                for i, prob in zip(valid_indices, probabilities):
                    results[i] = np.argmax(prob) if threshold_test(prob, threshold) else -1
            else:
                predictions = model.predict(valid_sequences)
                for i, pred in zip(valid_indices, predictions):
                    results[i] = pred
        else:
            if ovr in [0, 1, 2]:
                binary_classifier = model.estimators_[ovr]
                predictions = binary_classifier.predict(valid_sequences)
                for i, pred in zip(valid_indices, predictions):
                    results[i] = int(pred)
    
    return results

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

def get_threshold(precision, threshold_file, retrieve='threshold'):
    
    # Load the model if threshold_file is a string
    if isinstance(threshold_file, str):
        threshold_file = load_model(threshold_file)
        
    # Extract thresholds for each model
    model_thresholds = {
        name: [
            threshold_file[name][top][close(precision, threshold_file[name][top])][retrieve]
            for top in range(3)
        ]
        for name in threshold_file.keys()
    }
    
    return model_thresholds

def load_fasta(filepath):
    """
    Read sequences and their identifiers from a FASTA file.

    Args:
    - filepath (str): Path to the FASTA file.

    Returns:
    - headers (list): A list of sequence identifiers from the FASTA file.
    - sequences (list): A list of sequences from the FASTA file.
    """
    headers = []
    sequences = []
    with open(filepath, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            headers.append(record.id)
            sequences.append(str(record.seq))
    return headers, sequences
    
def read_fasta(data):
    """
    Parses the input data to extract sequences.
    If the data is in FASTA format, it extracts the sequence names and sequences.
    Otherwise, it treats the data as a simple sequence format.
    """
    sequences = []
    sequence_names = []
    temp_count = 0

    if ">" in data:  # Check if data is in FASTA format
        lines = data.strip().split("\n")
        for line in lines:
            if line.startswith(">"):
                sequence_names.append(line[1:])
                sequences.append("")
            else: 
                if not sequences:
                    temp_name = f"unnamed_sequence{temp_count}"
                    sequence_names.append(temp_name)
                    sequences.append("")
                    temp_count += 1
                sequences[-1] += line
                    
    else:
        sequences = data.strip().split("\n")

    return sequence_names, sequences

#%% GUI functions

def update_recall_and_threshold(precision, model_name):
    recall_list = [round(x, 2) for x in get_threshold(precision, 
                                                      os.path.join(application_path, '(precision)optimized_threshold.pkl'), 
                                                      retrieve='recall')[model_name]]
    threshold_list = [round(x, 2) for x in get_threshold(precision, 
                                                         os.path.join(application_path, '(precision)optimized_threshold.pkl'), 
                                                         retrieve='threshold')[model_name]]
    window['-RECALL-P-'].update(recall_list[0])
    window['-RECALL-AP-'].update(recall_list[1])
    window['-RECALL-H-'].update(recall_list[2])
            
    return recall_list, threshold_list

#%%
sg.LOOK_AND_FEEL_TABLE['CustomTheme'] = {'BACKGROUND': '#FFFFFF',
                                         'TEXT': '#000000',
                                         'INPUT': '#F2EFE8',
                                         'TEXT_INPUT': '#000000',
                                         'SCROLL': '#c7e78b',
                                         'BUTTON': ('white', '#709053'),
                                         'PROGRESS': ('#01826B', '#D0D0D0'),
                                         'BORDER': 1, 'SLIDER_DEPTH': 0, 'PROGRESS_DEPTH': 0}

sg.LOOK_AND_FEEL_TABLE['P2'] = {'BACKGROUND': '#e6f0f3',
                                           'TEXT': 'black',
                                           'INPUT': '#ffffff',
                                           'TEXT_INPUT': 'black',
                                           'SCROLL': '#d0e0e3',
                                           'BUTTON': ('black', '#b0c4c9'),
                                           'PROGRESS': ('#01826B', '#D0D0D0'),
                                           'BORDER': 1,
                                           'SLIDER_DEPTH': 0,
                                           'PROGRESS_DEPTH': 0}

sg.theme('P2')  

sg.set_options(font=('Helvetica', 16))

params_column = [
    [sg.Text("Model"), sg.Combo(['RandomForest (default)', 'CatBoostClassifier', 'ExtraTreesClassifier',
                                 'GradientBoostingClassifier', 'LightGBMClassifier', 'XGBClassifier'
                                 ], key='-MODEL-', default_value='RandomForest (default)', enable_events=True)],
    
    [sg.Text("Precision"), 
     sg.Spin(['{:.2f}'.format(i/100) for i in range(70, 96)], initial_value='default', key='-PRECISION-', size=(15, 1), enable_events=True)],
    [sg.Text("Recall (4+0)"), sg.Text('default', size = (5,1), key='-RECALL-P-')],
    [sg.Text("Recall (2+2)"), sg.Text('default', size = (5,1), key='-RECALL-AP-')],
    [sg.Text("Recall (3+1)"), sg.Text('default', size = (5,1), key='-RECALL-H-')],
]


params_frame = sg.Frame('', params_column, key='-PARAMS-', visible=False)

layout = [
    [sg.Text("Enter DNA sequence or upload .fasta file"),sg.Text('', size=(15,1), justification='right'), sg.Button('Show Example', key='-SHOW-EXAMPLE-')],
    [
        sg.Column([
            [sg.Multiline(key='-DNA-', size=(60, 5), expand_x=True, expand_y=True)],   # existing Multiline
            [sg.FileBrowse('Browse', key='-FILE-', initial_folder=os.getcwd())],
            [sg.Button('Submit', key='-SUBMIT-'), 
             sg.Button('Show Parameters', key='-SHOW_PARAMS-'),
             sg.Button('', image_filename=os.path.join(application_path, 'question_mark.png'), key='-QUESTION-MARK-', border_width=2,
                       pad=None, tooltip='Click for help', button_color=('black', '#e6f0f3'), visible=False),
             sg.Multiline(key='-RHS-SPACE-', size=(45, 5), expand_x=True, expand_y=True, visible=False,  font=("Helvetica", 11))
             ],
            [sg.pin(params_frame)],
            [sg.Button('Export results to csv', key='-EXPORT_CSV-', disabled=True), sg.Button('Reset', key='-RESET-')],
            [sg.Multiline(size=(60,9), key='-OUTPUT-', disabled=True, visible=False)]
        ])
    ]
]

window = sg.Window("G4 Predictor", layout, resizable=True, keep_on_top=True)

threshold_list = ''
model = load_model(os.path.join(application_path, 'RandomForest (default).pkl'))

# Event loop
while True:
    event, values = window.read()
        
    if event == sg.WINDOW_CLOSED or event == 'Quit':
        break
    
    if event == '-FILE-':
        if values['-FILE-']:
            if is_valid_fasta(values['-FILE-']):
                headers, sequences = load_fasta(values['-FILE-'])
                sequences = [seq.upper() for seq in sequences]
                combined = [f">{header}\n{seq}" for header, seq in zip(headers, sequences)]
                window['-DNA-'].update('\n'.join(combined))
            else:
                window['-OUTPUT-'].update('fasta file contains invalid characters (non A, T, C, G, N)', visible=True)

    if event == '-MODEL-':
        if values['-PRECISION-'] != 'default':
            precision = float(values['-PRECISION-'])
            recall_list, threshold_list = update_recall_and_threshold(precision, model_name=values['-MODEL-'])
            
        model = load_model(os.path.join(application_path, f"{values['-MODEL-']}.pkl"))
            
    if event == '-PRECISION-':
        precision = float(values['-PRECISION-'])
        recall_list, threshold_list = update_recall_and_threshold(precision, model_name=values['-MODEL-'])
    
    if event == '-RESET-':
        window['-DNA-'].update('')  
        window['-OUTPUT-'].update('', visible=False)  
        window['-RHS-SPACE-'].update('', visible=False)  
        window['-EXPORT_CSV-'].update(disabled=True)
        window['-MODEL-'].update(value='RandomForest (default)') 
        window['-PRECISION-'].update(value='default')  
        window['-RECALL-P-'].update(value='default')
        window['-RECALL-AP-'].update(value='default')
        window['-RECALL-H-'].update(value='default')
        

    if event == '-SHOW-EXAMPLE-':
        headers, sequences = load_fasta(os.path.join(application_path, 'sample.fasta'))
        combined = [f">{header}\n{seq}" for header, seq in zip(headers, sequences)]
        window['-DNA-'].update('\n'.join(combined))
    
    elif event == '-SHOW_PARAMS-':
        window['-PARAMS-'].update(visible=not window['-PARAMS-'].visible)
        btn_text = 'Hide Parameters' if window['-SHOW_PARAMS-'].get_text() == 'Show Parameters' else 'Show Parameters'
        window['-SHOW_PARAMS-'].update(btn_text)
        
        window['-QUESTION-MARK-'].update(visible=True) if window['-PARAMS-'].visible else window['-QUESTION-MARK-'].update(visible=False)
        
    if event == '-QUESTION-MARK-':
        with open(os.path.join(application_path, 'param_help.txt'), 'r') as file: ##
            text_content = file.read()
    
        current_content = window['-RHS-SPACE-'].get()
        if current_content == text_content:
            window['-RHS-SPACE-'].update(visible=False, value='')
        else:
            window['-RHS-SPACE-'].update(visible=True, value=text_content)

    if event == '-SUBMIT-':
        if not bool(values['-DNA-']):
            window['-OUTPUT-'].update('Please enter DNA sequence', visible=True)
        else:
            sequence_names, input_sequences = read_fasta(values['-DNA-'])
            original_sequences = [x.strip().upper() for x in input_sequences if x.strip() != '']
            
            g4_propensity = [G4Hunter(x) for x in original_sequences]
            sequences = seqs_converter(original_sequences)
            results = predict_G4(sequences, model=model, threshold=threshold_list)
            
            if values['-MODEL-'] == 'CatBoostClassifier':
                output = [f"{x} is {y}" for x, y in zip(original_sequences, 
                                                         [top_dict(inverted=True)[x.item()] for x in results])]
            else:
                output = [f"{x} is {y}" for x, y in zip(original_sequences, 
                                                        [top_dict(inverted=True)[x] for x in results])]
            
            if any('Invalid' in res for res in output):
                output.insert(0, "One or more sequences are invalid - only letters A, T, C, G, N are allowed.\n")
                
            window['-OUTPUT-'].update('\n'.join(output), visible=True)
            
            window['-EXPORT_CSV-'].update(disabled=False)
            
            warnings = []

            if any(x < 1.2 for x in g4_propensity):
                output_propensity = [f"{seq} may not form G4 (Invalid sequence)" if not validate_sequence(seq) else f"{seq} may not form G4" 
                                     for seq, prop in zip(original_sequences, g4_propensity) if prop < 1.2]
                warnings.append("WARNING. The following has low G4Hunter score:\n" + '\n'.join(output_propensity))
            
            if any(len(seq) <= 12 and validate_sequence(seq) for seq in original_sequences):
                output_length = [f"{seq} is too short" for seq in original_sequences if len(seq) <= 12 and validate_sequence(seq)]
                warnings.append("WARNING. The following sequences are too short:\n" + '\n'.join(output_length))
            
            if warnings:
                window['-RHS-SPACE-'].update('\n\n'.join(warnings), visible=True)
                    
    elif event == '-EXPORT_CSV-':
        df = pd.DataFrame([x.rsplit(' is ') for x in output])
        df.columns = ['Sequence', 'Topology']
        df['Comment'] = ['Sequence length too short' if len(x) <= 12 else 'G4Hunter Score too low' if G4Hunter(x) < 1.2 else '' for x in original_sequences]
        
        df['Model'] = [values['-MODEL-']] + ['' for x in range(len(original_sequences) - 1)]
        df['Precision'] = [values['-PRECISION-']] + ['' for x in range(len(original_sequences) - 1)]
        df['Recall_4+1'] = [window['-RECALL-P-'].DisplayText] + ['' for x in range(len(original_sequences) - 1)]
        df['Recall_2+2'] = [window['-RECALL-AP-'].DisplayText] + ['' for x in range(len(original_sequences) - 1)]
        df['Recall_3+1'] = [window['-RECALL-AP-'].DisplayText] + ['' for x in range(len(original_sequences) - 1)]
        if bool(sequence_names):
            df.insert(0, 'Name', sequence_names + ['' for x in range(len(original_sequences) - len(sequence_names))])

        save_path = sg.popup_get_file('Save results to CSV', save_as=True, default_path="g4_topology_prediction.csv", no_window=True, file_types=(("CSV Files", "*.csv"),))
        
        if save_path:
            df.to_csv(save_path, index=False)
        
    else:
        continue

# Close the window
window.close()


import random
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

AA_MASSES = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
    'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333
}

def cleave_sequence(sequence, num_cleaves):
    cleavage_sites = sorted(random.sample(range(1, len(sequence)), num_cleaves))
    annotated_cleavages = [list(range(i + 1, j + 1)) for i, j in zip([0] + cleavage_sites, cleavage_sites + [len(sequence)])]
    cleaved_sequences = [sequence[i:j] for i, j in zip([0] + cleavage_sites, cleavage_sites + [None])]
    return annotated_cleavages, cleaved_sequences

def get_ion_type_and_label(fragment, sequence_length):
    if fragment[0] == 1:
        return 'b', f'b{len(fragment)}'
    elif fragment[-1] == sequence_length:
        return 'y', f'y{len(fragment)}'
    else:
        return 'i', f'I{fragment[0]}-{fragment[-1]}'

def compute_mass(fragment_seq, ion_type):
    mass = sum(AA_MASSES[aa] for aa in fragment_seq)
    mass += 1.00727
    if ion_type == 'y':
        mass += 18.01560 - 1.00727
    return round(mass, 5)

def process_simulation_with_progress(args):
    num_simulations, num_cleaves, sequence = args
    sequence_length = len(sequence)
    ion_info = {'b': {}, 'y': {}, 'i': {}}

    for sim in tqdm(range(num_simulations), desc="Simulation"):
        annotated_cleavages, cleaved_sequences = cleave_sequence(sequence, num_cleaves)
        for fragment, fragment_seq in zip(annotated_cleavages, cleaved_sequences):
            ion_type, label = get_ion_type_and_label(fragment, sequence_length)
            mass = compute_mass(fragment_seq, ion_type)
            seq_str = ','.join(str(x) for x in fragment)
            if ion_type != 'i':
                key = (label, fragment_seq, mass)
                ion_info[ion_type][key] = ion_info[ion_type].get(key, 0) + 1
            else:
                key = mass
                if key not in ion_info['i']:
                    ion_info['i'][key] = {'Internal Fragment Label': label, 'Internal Fragment Sequence': fragment_seq, 'Internal Fragment Mass': mass, 'Internal Fragment Count': 0}
                else:
                    ion_info['i'][key]['Internal Fragment Label'] += ',' + label
                    ion_info['i'][key]['Internal Fragment Sequence'] += ',' + fragment_seq
                ion_info['i'][key]['Internal Fragment Count'] += 1

    return ion_info

def main(num_simulations, num_cleaves, sequence, run_index):
    num_cpus = cpu_count()
    args_list = [(num_simulations // num_cpus, num_cleaves, sequence)] * num_cpus

    with Pool(num_cpus) as pool:
        results = list(pool.imap(process_simulation_with_progress, args_list))

    combined_results = {'b': {}, 'y': {}, 'i': {}}
    for result in results:
        for ion_type in ('b', 'y', 'i'):
            for key, value in result[ion_type].items():
                if ion_type == 'i':
                    if key not in combined_results[ion_type]:
                        combined_results[ion_type][key] = value
                    else:
                        combined_results[ion_type][key]['Internal Fragment Label'] += ',' + value['Internal Fragment Label']
                        combined_results[ion_type][key]['Internal Fragment Sequence'] += ',' + value['Internal Fragment Sequence']
                        combined_results[ion_type][key]['Internal Fragment Count'] += value['Internal Fragment Count']
                else:
                    combined_results[ion_type][key] = combined_results[ion_type].get(key, 0) + value

    b_ions_df = pd.DataFrame([(key[0], key[1], key[2], count) for key, count in combined_results['b'].items()], columns=['b-ion Label', 'b-ion Sequence', 'b-ion Mass', 'b-ion Count']).sort_values(by='b-ion Mass')
    y_ions_df = pd.DataFrame([(key[0], key[1], key[2], count) for key, count in combined_results['y'].items()], columns=['y-ion Label', 'y-ion Sequence', 'y-ion Mass', 'y-ion Count']).sort_values(by='y-ion Mass')
    
    if combined_results['i']:
        internal_fragments_df = pd.DataFrame(list(combined_results['i'].values())).sort_values(by='Internal Fragment Mass')
    else:
        columns = ['Internal Fragment Label', 'Internal Fragment Sequence', 'Internal Fragment Mass', 'Internal Fragment Count']
        internal_fragments_df = pd.DataFrame(columns=columns)
    
    with pd.ExcelWriter('Observed_ion_frequencies.xlsx', engine='openpyxl') as writer:
        b_ions_df.to_excel(writer, sheet_name='b-ions', index=False)
        y_ions_df.to_excel(writer, sheet_name='y-ions', index=False)
        internal_fragments_df.to_excel(writer, sheet_name='Internal Fragments', index=False)
    filename = f'Observed_ion_frequencies_{run_index}.xlsx'
    with pd.ExcelWriter(filename, engine='openpyxl') as writer:
        b_ions_df.to_excel(writer, sheet_name='b-ions', index=False)
        y_ions_df.to_excel(writer, sheet_name='y-ions', index=False)
        internal_fragments_df.to_excel(writer, sheet_name='Internal Fragments', index=False)

    print(f"Results saved to '{filename}'")

if __name__ == "__main__":
    sequence = input("Enter the peptide sequence: ").upper()
    num_cleaves = int(input(f"Enter the number of cleavages you'd like to perform (max {len(sequence) - 1}): "))
    num_simulations = int(input(f"Enter the number of simulations you'd like to perform for each run: "))
    num_runs =  int(input(f"Enter the number of runs you'd like to perform : "))
    
    for run_index in range(1, num_runs+1): # Running 20 times
        print(f"Running simulation {run_index} of {num_runs}")
        main(num_simulations, num_cleaves, sequence, run_index)

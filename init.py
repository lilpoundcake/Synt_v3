import argparse
import math
import os
import numpy as np
import pandas as pd
import primer3
import plotly.io
import plotly.express as px

parser = argparse.ArgumentParser()

parser.add_argument('--name', '-n', default='synt_gene', help='Enter your project name')
parser.add_argument('--temp', '-t', type=int, default=65, help='Enter preferable temperature of annealing in PCR')
parser.add_argument('--primer_lenght', '-l', type=int, default=50, help='Select oligos length. Optimal length is 50-60 nt')
parser.add_argument('--expression_system', '-e', default='cho', choices=['cho', 'ecoli', 'sf9', 'human'], help='Select expression system')
parser.add_argument('--sequence_type', '-p', default='protein', choices=['protein', 'nucleotide'], help='Select type of input sequence')

# i'll add key for nogaps mode
# and codon usage distribution analysis tool

parser.add_argument('--codon_usage_threshold', '-u', default=20, type=float)
parser.add_argument('--tbio', default=True) 

parser.add_argument('seq', type=str)

args = parser.parse_args()

project_name = args.name
reaction_temp = args.temp
primer_lenght = args.primer_lenght
expression_system = args.expression_system
sequence_type = args.sequence_type
codon_usage_threshold = args.codon_usage_threshold

seq = args.seq



def split_seq(seq, fragment_num):
    fragment_list = []
    fragment_len = math.ceil(len(seq) / fragment_num)
    for i in range(fragment_num):
        fragment_list.append(seq[i*fragment_len:i*fragment_len + fragment_len])
    return fragment_list


if args.sequence_type == 'protein':
    len_check = 166
else:
    len_check = 500


fragment_list = []

if len(seq) > (len_check * 1.5):
    fragment_num = math.ceil(len(seq) / len_check)
    fragment_list = split_seq(seq, fragment_num)
else:
    fragment_list.append(seq)


os.mkdir(project_name)

def print_head(input_file, project_name, reaction_temp, primer_lenght, expression_system, codon_usage_threshold):
    print(f'melting low {reaction_temp}\nlength low {primer_lenght}\ntbio', file=input_file)
    print(f'frequency threshold {codon_usage_threshold}', file=input_file)
    print(f'logfile {project_name}_{i}.txt\n', file=input_file)

    if expression_system == 'cho' or expression_system == 'sf9':
        with open(f'codon_usage_tables/{expression_system}.txt', 'r') as codon_usage_table:
            print('codon', file=input_file)
            for line in codon_usage_table:
                print(f'{line.strip()}', file=input_file)
    elif expression_system == 'ecoli':
        print('codon E. coli\n', file=input_file)
    elif expression_system == 'human':
        print('codon H. sapiens\n', file=input_file)


for i, fragment in enumerate(fragment_list):
    with open(f'{project_name}/DNAWorks_{i}.inp', 'w') as input_file:
        print_head(input_file, project_name, reaction_temp, primer_lenght, 
                   expression_system, codon_usage_threshold)

        print(f'{sequence_type}', file=input_file)
        for j in range(math.ceil(len(fragment)/100)):
            print(f'{fragment[j*100:j*100+100]}', file=input_file)
        print('//\n', file=input_file)

for i, _ in enumerate(fragment_list):
    os.system(f'./DNAWorks/dnaworks {project_name}/DNAWorks_{i}.inp')

os.system(f'mv {project_name}_* {project_name}')



primer_dict = {}
primers_dict = {}
seq_dict = {}
seq = ''
primer_num = 0

for logfile_num, _ in enumerate(fragment_list):
    with open(f'{project_name}/{project_name}_{logfile_num}.txt', 'r') as log_file:
        for line in log_file:
            if "oligonucleotides need to be synthesized" in line:
                log_file.readline()
                line = log_file.readline()
                with open(f'{project_name}/{project_name}_{logfile_num}_primers.fa', 'w') as primers_fasta:
                    while " \n" != line:
                        primer_name = f'{project_name}_{primer_num}'
                        primer_seq = f'{line.split()[1]}'
                        primer_dict[primer_name] = primer_seq
                        print(f'>{primer_name}\n{primer_seq}', file=primers_fasta)
                        line = log_file.readline()
                        primer_num += 1
                primers_dict[logfile_num] = primer_dict
                primer_dict = {}

            elif "The DNA sequence" in line:
                log_file.readline()
                line = log_file.readline()
                with open(f'{project_name}/{project_name}_fragment_{logfile_num}_seq.fa', 'w') as seq_fasta:
                    while "---" not in line:
                        seq += line.split()[1]
                        line = log_file.readline()
                    seq_dict[logfile_num] = seq
                    print(f'>{project_name}_fragment_{logfile_num}\n{seq}', file=seq_fasta)
                    seq = ''

df1 = pd.DataFrame(primers_dict)

df2 = df1.reset_index().rename(columns={'index' : 'primer_name'}).melt(
    id_vars=['primer_name'], var_name='fragment', value_name='primer_seq'
    ).dropna(ignore_index=True)

df2['tm'] = df2.apply( lambda x : primer3.calc_tm(x.primer_seq), axis=1)
df2['hairpin_tm'] = df2.apply( lambda x : primer3.calc_hairpin_tm(x.primer_seq), axis=1)
df2['homodimer_tm'] = df2.apply( lambda x : primer3.calc_homodimer_tm(x.primer_seq), axis=1)
df2['gc_content'] = df2.apply( lambda x : (
    (x.primer_seq).count('G') + (x.primer_seq).count('C')
    ) / len(x.primer_seq) * 100, 
    axis=1)

df2.to_csv(f'{project_name}/{project_name}_primer_analysis.csv')
df2.query(f'hairpin_tm > {reaction_temp}').to_csv(f'{project_name}/{project_name}_high_hairpin_tm.csv')
df2.query(f'homodimer_tm > {reaction_temp}').to_csv(f'{project_name}/{project_name}_high_homodimer_tm.csv')



for fragment_num, _ in enumerate(fragment_list):
    df3 = df2.query(f'fragment == {fragment_num}')[['primer_name', 'primer_seq']]

    hetero_tm = {}
    high_heterodimer_tm = []
    high_temp_list = []

    for i in df3.primer_name:
        for j in df3.primer_name:
            hetero_tm[i] = hetero_tm.get(i, {})
            hetero_tm[i][j] = primer3.calc_heterodimer_tm(
                df3.query(f'primer_name == "{i}"').iat[0,1], 
                df3.query(f'primer_name == "{j}"').iat[0,1]
                )
            if hetero_tm[i][j] >= reaction_temp:
                if [j, i] in high_temp_list:
                    continue
                else:
                    high_heterodimer_tm.append([i, j, hetero_tm[i][j]])
                    high_temp_list.append([i, j])

    df4 = pd.DataFrame(hetero_tm)

    fig = px.imshow(pd.DataFrame(hetero_tm), color_continuous_scale='Viridis', text_auto=True, aspect="auto",
                    title=f'Fragment_{fragment_num}', labels=dict(x="Primer_1", y="Primer_2", color='Tm'))
    
    plotly.io.write_html(fig, f'{project_name}/{project_name}_fragment_{fragment_num}_heterodimer_heatmap.html')

    df4.to_csv(f'{project_name}/{project_name}_fragment_{fragment_num}_heterodimer_tm.csv')

    pd.DataFrame(
        high_heterodimer_tm, 
        columns=['primer_1', 'primer_2', 'tm']
        ).to_csv(f'{project_name}/{project_name}_fragment_{fragment_num}_high_heterodimer_tm.csv')


for fragment_num, seq_name in enumerate(seq_dict):
    fragment_seq = seq_dict[seq_name]
    gc_list = []
    gc_min = float(100)
    gc_max = float()
    gc_mean = float()

    for i in range(len(fragment_seq)-29):
        frame = fragment_seq[i:i+30]
        gc_content = (frame.count('C') + frame.count('G')) / len(frame) * 100
        if gc_content > gc_max:
            gc_max = gc_content
        elif gc_content < gc_min:
            gc_min = gc_content
        gc_mean += gc_content
        gc_list.append(gc_content)
    gc_mean = gc_mean / len(gc_list)

    df5 = pd.DataFrame({'gc' : gc_list,
             'gc_max' : gc_max, 
             'gc_min' : gc_min,
             'gc_mean' : gc_mean})
    
    fig = px.line(df5, labels=dict(index='Rollmean_30nt', value='GC_content'), 
              color_discrete_sequence=px.colors.qualitative.Pastel1,
              title=f'Fragment {fragment_num}. Rollmean GC-content. Frame 30 nt')
    
    plotly.io.write_html(fig, f'{project_name}/{project_name}_fragment_{fragment_num}_rollmean_gc_content.html')

full_seq = ''

with open(f'{project_name}/{project_name}_synthetic_gene.fa', 'w') as gene:
    print(f'>{project_name}_synthetic_gene', file=gene)
    for seq in seq_dict:
        full_seq += seq_dict[seq]
    for i in range(math.ceil(len(full_seq)/100)):
        print(full_seq[i*100:i*100+100], file=gene)

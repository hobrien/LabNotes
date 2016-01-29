from os import listdir, path

folder = 'LingulodiniumBlast'
def GetGenes(file): #Extract gene name from file name
    return file.split('_')[0]

def task_extract_top_hit():
    for file in listdir(folder):
        if '_top_hit' in file:
            continue
        gene = GetGenes(file)
        yield {
            'name': gene,
            'file_dep': [ path.join(folder,file) ],
            'targets': [ path.join(folder, '{0}_top_hit.bl'.format(gene)) ],
            'actions': ['head -1 %(dependencies)s > %(targets)s']
        }

def task_extract_top_hit_sequence():
    for file in listdir(folder):
        if '_top_hit.bl' not in file:
            continue
        gene = GetGenes(file)
        yield {
            'name': gene,
            'file_dep': [ path.join(folder,file) ],
            'targets': [ path.join(folder, '{0}_top_hit.fa'.format(gene)) ],
            'actions': ["ParseBlast.py -p tblastn -t --outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore' %(dependencies)s Lingulodinium.fa > %(targets)s"]
        }

def task_blast_top_hit_sequence():
    for file in listdir(folder):
        if '_top_hit.fa' not in file:
            continue
        gene = GetGenes(file)
        yield {
            'name': gene,
            'file_dep': [ path.join(folder,file) ],
            'targets': [ path.join('LingulodiniumSwissProt', '{0}_Lingulodinium.fa'.format(gene)) ],
            'actions': ["blastp -remote -query %(dependencies)s -db swissprot -max_target_seqs 1  -outfmt '6 qseqid sacc stitle evalue' | head -1  > %(targets)s"]
        }

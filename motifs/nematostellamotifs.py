from Bio import SeqIO
import re
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm
import sys
import argparse
import os

import warnings
warnings.filterwarnings('ignore')

"""
Version 2 of the script I'm writing to find GO terms associated with motifs in the Nematostella genome.
Changes here include the addition of progress bars, usage of dataframes rather than dictionaries
for later analysis of specific elements, and the removal of duplicates done with these dataframes. 

Currently I'm struggling with my approach to this, but I think I have the collection of artifacts down.
What I want to do next is optimize the collection of GO terms from said artifacts, for which I think the
best strategy is to do (for each row in artifacts): extract geneID -> get corresponding GO terms with loc (not iloc) -> add
to dictionary -> rinse and repeat.
"""

#read nematostella genome to be accessible to all functions
genome = []
with open("../Nvec200.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        genome.append(record)

#read nematostella GTF, also accessible to all functions
colnames = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf = pd.read_table("../../NV2g.20240221.gtf", names = colnames)

colnames_go = ["id", "type", "term_id", "term_name", "description", "source"]
goterms = pd.read_table("../../NV2g.20240221.GO.tsv", names=colnames_go)

chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
        "chrUn1", "chrUn2", "chrUn3", "chrUn4", "chrUn5", "chrUn6", "chrUn7", "chrUn8", "chrUn9", "chrUn10", "chrUn11", "chrUn12", "chrUn13", "chrUn14", "chrUn15"]

def find_motif_positions(motif, ref):
    """
    Given a motif string and a reference SeqRecord object from BioPython, find all the starting positions of
    that motif in the SeqRecord. Ignores partial or imperfect alignments. Returns a list of all starting positions.
    - motif: the actual string that you want to find
    - ref: a SeqRecord object that may or may not contain the motif.
    """
    #assumes ref is a SeqRecord object; we convert it to a string for using RegularExpressions
    ref_string = str(ref.seq)
    return [m.start() for m in re.finditer(motif, ref_string)]

def find_closest_upstream(pos, chr, within = 10000):
    """
    Given a chromosomal position and a chromosome name, this function will find the closest features downstream of the position
    within a certain radius; that is, having a higher positional start value than the position parameter. 
    Objects are returned as a pandas dataframe.
    - pos: integer position for which you want to find downstream object(s).
    - chr: name of the chromosome on which you want to find object(s).
    - within (default 10000): integer number for the distance from start site you wish to examine.
    """
    df_filt = gtf[gtf["chr"] == chr]
    df_filt = df_filt[df_filt["start"] <= pos]
    df_filt = df_filt[df_filt["start"] >= (pos - within)]

    #the GTF is NOT sorted by default, so we'll sort the GTF subset
    df_filt_sort = df_filt.sort_values(by=["start", "end"], ascending=False)
    return df_filt_sort

def find_closest_downstream(pos, chr, within = 10000):
    """
    Given a chromosomal position and a chromosome name, this function will find the closest features downstream of the position
    within a certain radius; that is, having a higher positional start value than the position parameter. 
    Objects are returned as a pandas dataframe.
    - pos: integer position for which you want to find downstream object(s).
    - chr: name of the chromosome on which you want to find object(s).
    - within (default 10000): integer number for the distance from start site you wish to examine.
    """
    df_filt = gtf[gtf["chr"] == chr]
    df_filt = df_filt[df_filt["start"] >= pos]
    df_filt = df_filt[df_filt["start"] <= (pos + within)]

    #the GTF is NOT sorted by default, so we'll sort the GTF subset
    df_filt_sort = df_filt.sort_values(by=["start", "end"])
    return df_filt_sort

def remove_id_term_dups(df):
    """
    Given a dataframe of elements, this function will find the unique id/term pairs and preserve only the first. 
    In the context of the GO TSV, instances of duplicate id/term pairs often exist primarily due to different sources
    (InterPro/OMA) and are functionally the same otherwise. In the context of collecting terms proximal to motif positions,
    instances of duplicate id/term pairs are double-counted and should have all but one removed.

    Later edit: This can easily be done with df.drop_duplicates() and so is not really needed.  
    """
    return df.drop_duplicates(subset=["id", "term_id"])

def find_goterms_onepos(pos, chr, within=10000):
    """
    Finds the GO terms linked to one specific position on a character. Meant to be used as a helper function.
    """
    objs_upstream = find_closest_upstream(pos, chr, within=within)
    objs_downstream = find_closest_downstream(pos, chr, within=within)

    all_objs = pd.concat([objs_downstream, objs_upstream])

    term_ids = {}
    term_names = {}
    for i in range(len(all_objs)):
        current_id = all_objs.iloc[i]["attribute"][42:-2]
        current_row = goterms[goterms["id"] == current_id]

        #if there is none: keep going
        if len(current_row) <= 0:
            continue
        else:
            current_row = current_row.iloc[0]

        if str(current_row["term_id"]) in term_ids:
            term_ids[str(current_row["term_id"])] += 1
        else:
            term_ids[str(current_row["term_id"])] = 1

        if str(current_row["term_name"]) in term_names:
            term_names[str(current_row["term_name"])] += 1
        else:
            term_names[str(current_row["term_name"])] = 1
    return term_ids, term_names

def find_goterms_onechr(positions, chr, within=10000):
    """ 
    Finds the GO terms associated with all the motif positions on a single chromosome. Also can meant to be used as a helper function.
    IMPORTANT: This function may count duplicates. 
    TODO: fix this function and previous helper ones so that duplicates are not counted. 
    """
    term_ids = {}
    term_names = {}
    for p in positions:
        #print(p)
        current_ids, current_names = find_goterms_onepos(p, chr, within=within)
        for key in current_ids:
            if key in term_ids:
                term_ids[key] += current_ids[key]
            else:
                term_ids[key] = current_ids[key]

        for name in current_names:
            if name in term_names:
                term_names[name] += current_names[name]
            else:
                term_names[name] = current_names[name]
    return term_ids, term_names

def find_artifacts_onepos(pos, chr, within=10000):
    """
    Finds the artifacts linked to one specific position on a character. Meant to be used as a helper function.
    """
    objs_upstream = find_closest_upstream(pos, chr, within=within)
    objs_downstream = find_closest_downstream(pos, chr, within=within)
    all_objs = pd.concat([objs_downstream, objs_upstream])

    # goterms_onepos = pd.DataFrame()

    # for i in range(len(all_objs)):
    #     current_id = all_objs.iloc[i]["attribute"][42:-2]
    #     current_row = goterms[goterms["id"] == current_id]

    #     #if there is none: keep going
    #     if len(current_row) <= 0:
    #         continue
    #     else:
    #         current_row = current_row.iloc[0]

    #     goterms_onepos = pd.concat([goterms_onepos, current_row])

    return all_objs

def find_artifacts_onechr(positions, chr, within=10000):
    """ 
    Finds the artifacts associated with all the motif positions on a single chromosome. Also can meant to be used as a helper function.
    IMPORTANT: This function on its own counts duplicates, which are removed with a different function later. 
    """
    terms = pd.DataFrame()
    for p in positions:
        current_terms = find_artifacts_onepos(p, chr, within=within)
        terms = pd.concat([terms, current_terms])

    return terms

def plot_goterms(counts, n=10, motif=""):
    """
    Takes a dataframe of counts, which can be obtained by as_dataframe(), and makes a horizontal barplot
    of the counts for the top n GO terms.
    - n (default 10): number of GO terms to display.
    """
    #does not assume counts are sorted; sorts them if so
    counts_sorted = counts.sort_values(by=["count"], ascending=False)
    counts_subframe = counts_sorted[0:n]

    plt.barh(counts_subframe['name'], counts_subframe["count"])
    plt.title('Enriched GO Term Counts')
    plt.ylabel('GO Term')
    plt.xlabel('Frequency')
    print("Saving term plot to " + motif + "_goterms.png")
    plt.savefig(motif+"_goterms.png")

def as_dataframe(counts):
    """
    Takes a dictionary of counts obtained from finding GO terms and returns a dataframe of the names/ids and respective counts.
    """
    name_list = list(counts)
    vals = [counts[n] for n in name_list]
    counts_df = pd.DataFrame({"name": name_list, "count": vals})
    return counts_df

def as_count_dict(df, col="term_name"):
    """
    Takes a dataframe of GO terms and finds the counts associated with each unique GO term, returned as a dictionary.
    Can also do this with term_id if you want, but as this function is used for plotting, we use the term names.
    """
    counts = {}
    for i in range(len(df)):
        key = df.iloc[i][col]
        if key in counts:
            counts[key] += 1
        else:
            counts[key] = 1
    return counts

def find_goterms(motif, within=10000):
    """
    Function that puts it all together to find all of the GO terms associated with a motif. Here, this takes the motif as a string,
    gets the positions along a chromosome, finds the GO terms from that chromosome, and repeats across the entire genome.
    - motif: string representation of the motif, for which exact matches are found.
    """

    assert len(chrs) == len(genome) #optional since we load these ourselves and I know it works

    artifacts = pd.DataFrame()

    print("Finding artifacts associated with motif " + motif)
    for i in tqdm(range(len(chrs))):
        assert genome[i].name == chrs[i] #important to ensure we keep positions consistent
        positions = find_motif_positions(motif, genome[i]) 

        current_terms = find_artifacts_onechr(positions, chrs[i], within=within)
        artifacts = pd.concat([artifacts, current_terms])

    print("Finding GO terms for artifacts for motif " + motif)
    for i in tqdm(range(len(chrs))):
        assert genome[i].name == chrs[i]
        positions = find_motif_positions(motif, genome[i]) 

        current_ids, current_names = find_goterms_onechr(positions, chrs[i], within=within)

        term_ids = {}
        term_names = {}
        for key in current_ids:
            if key in term_ids:
                term_ids[key] += current_ids[key]
            else:
                term_ids[key] = current_ids[key]

        for name in current_names:
            if name in term_names:
                term_names[name] += current_names[name]
            else:
                term_names[name] = current_names[name] 

    return artifacts, term_ids, term_names

def process_motif(motif, within=10000, n=10):
    """
    Full function to process a singular motif, from getting the full deduplicated dataframe of GO terms associated with it to
    making the plot of the top n terms. Both the dataframe and the plot are saved. 
    """

    artifacts, terms_ids, term_names = find_goterms(motif, within=within)
    #terms_df = remove_id_term_dups(terms_df)

    print("Saving artifacts to " + motif + "_artifacts.csv")
    artifacts.to_csv(motif + "_artifacts.csv")

    counts = as_count_dict(term_names)
    counts_df = as_dataframe(counts)

    plot_goterms(counts_df, n=n, motif=motif)

def parse_args():
    """
    Standard argument parsing to pass either a singular motif or a txt file containing motifs (one per line) to the 
    script via command line.
    """
    parser = argparse.ArgumentParser(description="Find GO terms associated with sequence motifs in the Nematostella vectensis genome.")
    parser.add_argument("--motif", required=False, help="Singular DNA sequence motif to be analyzed")
    parser.add_argument("--file", required=False, help="Path to text file containing multiple sequence motifs to be analyzed (one per line).")
    parser.add_argument("--within", required=False, help="Distance from motif to check for GO terms. Default, 10kpb.")
    parser.add_argument("--n", required=False, help="Number of top associated GO terms to plot. Default, 10.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    motif = ""
    if args.motif is not None:
        motif = args.motif

    motifs = []
    if args.file is not None:
        assert os.path.isfile(args.file)
        with open(args.file, newline="") as fr:
            for line in fr:
                motifs.append(line.strip())

    within = 10000
    if args.within is not None:
        within = int(args.within)

    n = 10
    if args.n is not None:
        n = args.n

    #case 1: we didn't find anything
    if len(motif) == 0 and len(motifs) == 0:
        print("No motifs passed, exiting program. Use -h or --help for more detailed instructions.")
        exit()

    #case 2: motif is of length > 0 AND motifs is of length > 0
    #in this case, we demand the user only do one
    if len(motif) > 0 and len(motifs) > 0:
        print("Please pass only one of --motif or --file to the program.")
        exit()

    #case 3: motif is of length > 0; at this point we will analyze a singular motif.
    if len(motif) > 0:
        process_motif(motif, within=within, n=n)

    #case 4: motifs is of length > 0; analyze all motifs in a loop.
    if len(motifs) > 0:
        for m in motifs:
            process_motif(m, within=within, n=n)

    else:
        print("If you've reached this point in the program, something went horribly. Congrats on messing this up beyond belief.")
        exit()

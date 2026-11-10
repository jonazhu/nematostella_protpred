from Bio import SeqIO
import re
import pandas as pd
import matplotlib.pyplot as plt

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

def find_goterms_onepos(pos, chr, within=10000):
    """
    Finds the GO terms linked to one specific position on a character. Meant to be used as a helper function.
    """
    objs_upstream = find_closest_upstream(pos, chr, within=within)
    objs_downstream = find_closest_downstream(pos, chr, within=within)

    all_objs = pd.concat([objs_downstream, objs_upstream])

    term_ids = []
    term_names = []
    for i in range(len(all_objs)):
        current_id = all_objs.iloc[i]["attribute"][42:-2]
        current_row = goterms[goterms["id"] == current_id]

        #if there is none: keep going
        if len(current_row) <= 0:
            continue
        else:
            current_row = current_row.iloc[0]

        term_ids.append(str(current_row["term_id"]))
        term_names.append(str(current_row["term_name"]))
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

def plot_goterms(counts, n=10):
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
    plt.show()

def as_dataframe(counts):
    """
    Takes a dictionary of counts obtained from finding GO terms and returns a dataframe of the names/ids and respective counts.
    """
    name_list = list(counts)
    vals = [counts[n] for n in name_list]
    counts_df = pd.DataFrame({"name": name_list, "count": vals})
    return counts_df

def find_goterms(motif, within=10000):
    """
    Function that puts it all together to find all of the GO terms associated with a motif. Here, this takes the motif as a string,
    gets the positions along a chromosome, finds the GO terms from that chromosome, and repeats across the entire genome.
    - motif: string representation of the motif, for which exact matches are found.
    """

    assert len(chrs) == len(genome)

    term_ids = {}
    term_names = {}

    for i in range(len(chrs)):
        assert genome[i].name == chrs[i]
        positions = find_motif_positions(motif, genome[i]) 

        current_ids, current_names = find_goterms_onechr(positions, chrs[i], within=within)
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
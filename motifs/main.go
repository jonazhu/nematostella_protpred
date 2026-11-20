package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"strings"
	"time"
)

var distance int
var numProcs int
var precision float64
var singularMotif string
var infile string
var genome string
var countfile string
var artifactfile string
var gotermfile string
var outfolder string

func init() {
	flag.IntVar(&distance, "w", 1000, "Distance in base pairs to check motifs for associated artifacts.")
	flag.IntVar(&numProcs, "numProcs", runtime.NumCPU(), "Number of processors to use for parallel processing of artifacts and GO terms. The default will be the number of processors you have, shown in the provided default.")
	flag.Float64Var(&precision, "p", 0.00001, "Float precision to compare floats.")
	flag.StringVar(&singularMotif, "m", "none", "Motif to look for in the genome; will be capitalized by the program automatically. If not specified, a text file containing one motif per line must be specified with -infile.")
	flag.StringVar(&infile, "f", "none", "Text file containing one motif per line to look for in the genome; motifs will be automatically capitalized. If not specified, a singular motif must be specified with -motif.")
	flag.StringVar(&genome, "g", "none", "Build of the genome to analyze for motifs: must be one of either Nvec200 or jaNemVect.")
	flag.StringVar(&countfile, "countFile", "counts.csv", "Base output filename for the counts of GO terms found.")
	flag.StringVar(&artifactfile, "artifactFile", "artifacts.csv", "Base output filename for the proximal artifacts found to motifs.")
	flag.StringVar(&gotermfile, "goTermFile", "assoc_goterms.csv", "Base output filename for the GO terms matched to proximal artifacts.")
	flag.StringVar(&outfolder, "o", ".", "Folder name for results to go into. Default is the current directory.")
}

func main() {
	flag.Parse()
	if (genome != "Nvec200") && (genome != "jaNemVect") {
		flag.Usage()
		os.Exit(1)
	}
	if singularMotif == "none" && infile == "none" {
		flag.Usage()
		os.Exit(1)
	}
	fmt.Println("Motif Analysis for Nematostella vectensis Genome")

	var background dataframe
	var artifacts dataframe

	//If Nvec200 is the specified genome, the following section is ran.
	//jaNemVect files are specified differently and so call different functions.
	if genome == "Nvec200" {
		gtf := ReadGTF("../data/NV2g.20240221.gtf")
		f := ReadFasta("../data/Nvec200.fasta")
		goterms := ReadGOTerms("../data/NV2g.20240221.GO.tsv")

		if !DoesFileExist("background_Nvec.csv") {
			fmt.Println("Background file not found. Making background count file.")
			bg_counts := goterms.CountGOTerms(3)
			background = CountsToDF(bg_counts)
			bg_ratios := background.FindTermRatios()
			background = background.AddNumericColumn(bg_ratios, "term_ratio")
			WriteCSV(background, "background_Nvec.csv")
		}

		background = ReadCSV("background_Nvec.csv")
		keys := GetKeys(f)

		motifs := make([]string, 1)
		if singularMotif != "none" {
			motifs[0] = strings.ToUpper(singularMotif)
		} else { //due to the usage statement we added earlier, this means we have a file of motifs
			motifs = ReadMotifsFromFile(infile)
		}

		for m := range motifs {
			currentMotif := motifs[m]
			fmt.Println("-------------------------------------")
			fmt.Println("Currently analyzing motif", currentMotif)

			start := time.Now()
			duration := time.Since(start)
			fmt.Println("Finding artifacts.")
			for k := range keys {
				//find the artifacts associated with each chromosome
				motif_pos := FindMotifPosOneChrom(currentMotif, f[keys[k]])
				artifacts_current := FindArtifactsOneChrom(motif_pos, gtf, keys[k][1:], distance, numProcs) //remove >
				artifacts = append(artifacts, artifacts_current...)

				duration = time.Since(start)
			}
			fmt.Println(len(artifacts), "artifacts found in", duration.Seconds(), "seconds")

			WriteCSV(artifacts, outfolder+"/"+currentMotif+"_"+artifactfile)

			start = time.Now()
			fmt.Println("Finding GO terms.")
			results := artifacts.FindGOTerms(goterms, numProcs, genome)
			duration = time.Since(start)
			fmt.Println(len(results), "GO terms found in", duration.Seconds(), "seconds")
			WriteCSV(results, outfolder+"/"+currentMotif+"_"+gotermfile)

			start = time.Now()
			fmt.Println("Counting GO terms.")
			counts := results.CountGOTerms(3)
			duration = time.Since(start)
			fmt.Println(len(counts), "GO terms counted in", duration.Seconds(), "seconds")
			counts_df := CountsToDF(counts)

			start = time.Now()
			fmt.Println("Processing GO terms.")
			ratios := counts_df.FindTermRatios()
			counts_df = counts_df.AddNumericColumn(ratios, "term_ratio")
			test_res := TestTermRatios(counts_df, background)
			counts_df = counts_df.AddNumericColumn(test_res[0], "z_score")
			counts_df = counts_df.AddNumericColumn(test_res[1], "p_val")
			padj := CorrectPVals(test_res[1], precision)
			counts_df = counts_df.AddNumericColumn(padj, "padj")
			duration = time.Since(start)
			fmt.Println(len(counts_df), "GO terms processed in", duration.Seconds(), "seconds")
			WriteCSV(counts_df, outfolder+"/"+currentMotif+"_"+countfile)
		}
	} else if genome == "jaNemVect" {
		gtf := ReadGTF("../data/genomic.gtf")
		f := ReadFasta("../data/GCF_932526225.1_jaNemVect1.1_genomic.fna")
		goterms := ReadOntology("../data/GCF_932526225.1_jaNemVect1.1_gene_ontology.gaf")

		term2name := ReadTerm2Name("../data/GOterm2name.csv")

		if !DoesFileExist("background_jaNemVect.csv") {
			fmt.Println("Background file not found. Making background count file.")
			bg_counts := goterms.CountGOTerms(4)
			background = CountsToDF(bg_counts)
			bg_ratios := background.FindTermRatios()
			background = background.AddNumericColumn(bg_ratios, "term_ratio")
			WriteCSV(background, "background_jaNemVect.csv")
		}

		background = ReadCSV("background_jaNemVect.csv")
		keys := GetKeys(f)

		motifs := make([]string, 1)
		if singularMotif != "none" {
			motifs[0] = strings.ToUpper(singularMotif)
		} else { //due to the usage statement we added earlier, this means we have a file of motifs
			motifs = ReadMotifsFromFile(infile)
		}

		for m := range motifs {
			currentMotif := motifs[m]
			fmt.Println("-------------------------------------")
			fmt.Println("Currently analyzing motif", currentMotif)

			start := time.Now()
			duration := time.Since(start)
			fmt.Println("Finding artifacts.")
			for k := range keys {
				//find the artifacts associated with each chromosome
				motif_pos := FindMotifPosOneChrom(currentMotif, strings.ToUpper(f[keys[k]]))
				artifacts_current := FindArtifactsOneChrom(motif_pos, gtf, keys[k][1:12], distance, numProcs) //remove > and trailing description
				artifacts = append(artifacts, artifacts_current...)

				duration = time.Since(start)
			}
			fmt.Println(len(artifacts), "artifacts found in", duration.Seconds(), "seconds")
			WriteCSV(artifacts, outfolder+"/"+currentMotif+"_"+artifactfile)

			start = time.Now()
			fmt.Println("Finding GO terms.")
			results := artifacts.FindGOTerms(goterms, numProcs, genome)
			duration = time.Since(start)
			fmt.Println(len(results), "GO terms found in", duration.Seconds(), "seconds")
			WriteCSV(results, outfolder+"/"+currentMotif+"_"+gotermfile)

			start = time.Now()
			fmt.Println("Counting GO terms.")
			counts := results.CountGOTerms(4)
			duration = time.Since(start)
			fmt.Println(len(counts), "GO terms counted in", duration.Seconds(), "seconds")
			counts_df := CountsToDF(counts)

			start = time.Now()
			fmt.Println("Processing GO terms.")
			ratios := counts_df.FindTermRatios()
			counts_df = counts_df.AddNumericColumn(ratios, "term_ratio")
			test_res := TestTermRatios(counts_df, background)
			counts_df = counts_df.AddNumericColumn(test_res[0], "z_score")
			counts_df = counts_df.AddNumericColumn(test_res[1], "p_val")
			padj := CorrectPVals(test_res[1], precision)
			counts_df = counts_df.AddNumericColumn(padj, "padj")

			//change the GO terms to the actual names in counts_df
			for i := 1; i < len(counts_df); i++ {
				current_term := counts_df[i][0]
				counts_df[i][0] = term2name[current_term]
			}

			duration = time.Since(start)
			fmt.Println(len(counts_df), "GO terms processed in", duration.Seconds(), "seconds")
			WriteCSV(counts_df, outfolder+"/"+currentMotif+"_"+countfile)
		}
	}
}

package main

import (
	"fmt"
	"time"
)

func main() {
	fmt.Println("Motif Analysis for Nematostella vectensis Genome")

	gtf := ReadGTF("../data/NV2g.20240221.gtf")
	f := ReadFasta("../data/Nvec200.fasta")
	goterms := ReadGOTerms("../data/NV2g.20240221.GO.tsv")

	var background dataframe

	if !DoesFileExist("background.csv") {
		fmt.Println("Background file not found. Making background count file.")
		bg_counts := goterms.CountGOTerms()
		background = CountsToDF(bg_counts)
		bg_ratios := background.FindTermRatios()
		background = background.AddNumericColumn(bg_ratios, "term_ratio")
		WriteCSV(background, "background.csv")
	}

	background = ReadCSV("background.csv")
	keys := GetKeys(f)

	var artifacts dataframe
	start := time.Now()
	fmt.Println("Finding artifacts.")
	for k := range keys {
		//find the artifacts associated with each chromosome
		motif_pos := FindMotifPosOneChrom("AATAAAAAAA", f[keys[k]])
		artifacts_current := FindArtifactsOneChrom(motif_pos, gtf, keys[k][1:], 1000, 4) //remove >
		artifacts = append(artifacts, artifacts_current...)
	}

	duration := time.Since(start)
	// Nanoseconds as int64
	fmt.Println(len(artifacts), "artifacts found in", duration.Seconds(), "seconds")

	WriteCSV(artifacts, "artifacts.csv")

	start = time.Now()

	fmt.Println("Finding GO terms.")
	results := artifacts.FindGOTerms(goterms, 4)
	duration = time.Since(start)
	fmt.Println(len(results), "GO terms found in", duration.Seconds(), "seconds")
	WriteCSV(results, "assoc_goterms.csv")

	start = time.Now()
	fmt.Println("Counting GO terms.")
	counts := results.CountGOTerms()
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
	padj := CorrectPVals(test_res[1], 0.00001)
	counts_df = counts_df.AddNumericColumn(padj, "padj")
	duration = time.Since(start)
	fmt.Println(len(counts_df), "GO terms processed in", duration.Seconds(), "seconds")
	WriteCSV(counts_df, "counts.csv")
}

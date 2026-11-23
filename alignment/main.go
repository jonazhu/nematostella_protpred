package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"runtime"
	"time"
)

var substitute string
var subLetter byte
var numProcs int
var fasta1 string
var fasta2 string
var subMatFile string
var outfasta string
var outpretty string
var outmatches string
var startGap int
var extendGap int
var descWidth int
var prettyWidth int

func init() {
	flag.StringVar(&substitute, "subU", "X", "Letter to substitute for U (selenocysteine), which appears in some protein sequences in jaNemVect1.1.")
	flag.IntVar(&numProcs, "numProcs", runtime.NumCPU(), "Number of processors to use for parallel processing of alignments. The default will be the number of processors you have, shown in the provided default.")
	flag.StringVar(&fasta1, "f1", "none", "Input FASTA file of proteins to be aligned.")
	flag.StringVar(&fasta2, "f2", "none", "Input FASTA file of proteins to which F1 proteins will be aligned.")
	flag.StringVar(&outfasta, "of", "alignments.fasta", "Output FASTA file of optimal alignments generated without scores.")
	flag.StringVar(&outpretty, "op", "alignments_pretty.txt", "Output file of optimal alignments neatly printed with scores and identity match proportions.")
	flag.StringVar(&outmatches, "om", "matches.csv", "Output file of best match pairings and the final alignment scores in CSV format.")
	flag.StringVar(&subMatFile, "sm", "blosum62.txt", "Text file containing the substitution matrix to be used in alignment. By default, the provided BLOSUM62 matrix is used, consistent with BLASTP.")
	flag.IntVar(&startGap, "sg", -12, "Penalty to be applied for starting a gap. MUST BE NEGATIVE. By default, the BLASTP penalty is used.")
	flag.IntVar(&extendGap, "eg", -1, "Penalty to be applied for extending a gap. MUST BE NEGATIVE. By default, the BLASTP penalty is used.")
	flag.IntVar(&descWidth, "dw", 15, "Number of characters to be used in the descriptions of sequences for pretty alignment file.")
	flag.IntVar(&prettyWidth, "pw", 150, "Number of alignment characters to be written to each line for the pretty alignment.")
}

func main() {
	flag.Parse()
	if fasta1 == "none" || fasta2 == "none" {
		flag.Usage()
		os.Exit(1)
	}
	fastas_Nvec200 := ReadFasta(fasta1)
	fastas_jaNemVect := ReadFasta(fasta2)
	keys_Nvec200 := GetKeys(fastas_Nvec200)
	keys_jaNamVect := GetKeys(fastas_jaNemVect)

	mat, mapping := ReadBlosum(subMatFile)
	mat, mapping = AddAlias(mat, mapping, 'U', substitute[0]) //by default we score selenocysteine as X, like BLAST does

	file, err := os.Create(outfasta) // truncates once, at the start of the run
	if err != nil {
		panic(err)
	}
	defer file.Close()

	w := bufio.NewWriter(file)
	defer w.Flush() // deferred last-in-first-out, so this runs before file.Close()

	prettyfile, err := os.Create(outpretty)
	if err != nil {
		panic(err)
	}
	defer prettyfile.Close()

	pw := bufio.NewWriter(prettyfile)
	defer pw.Flush()

	matchfile, err := os.Create(outmatches)
	if err != nil {
		panic(err)
	}
	defer matchfile.Close()

	mw := bufio.NewWriter(matchfile)
	defer mw.Flush()
	fmt.Fprintln(mw, "protein1,protein2,score")

	for k1 := range keys_Nvec200 {
		start := time.Now()
		fmt.Println("Aligning", keys_Nvec200[k1])
		p1 := fastas_Nvec200[keys_Nvec200[k1]]
		best_align := FindBestAlignment(p1, keys_Nvec200[k1], keys_jaNamVect, fastas_jaNemVect, mat, mapping, startGap, extendGap, numProcs)
		fmt.Println("Best match:", best_align.name, "with a score of", best_align.score)
		s1, s2, _ := Traceback(p1, fastas_jaNemVect[best_align.name], best_align.alignmentMatrix, startGap, extendGap)
		WriteAlignment(w, keys_Nvec200[k1], best_align.name, s1, s2)
		fmt.Fprintln(w, "")
		WritePrettyAlignment(pw, keys_Nvec200[k1][0:descWidth+1], best_align.name[0:descWidth+1], s1, s2, best_align.score, prettyWidth)
		fmt.Fprintln(pw, "")
		WriteMatch(mw, keys_Nvec200[k1], best_align.name, best_align.score)

		duration := time.Since(start)
		fmt.Println("Finished processing", keys_Nvec200[k1], "in", duration)
		fmt.Println("")
	}
	fmt.Println("All proteins aligned")
}

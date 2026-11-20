package main

import (
	"bufio"
	"errors"
	"fmt"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
)

/*
-------------------------------
DATATYPE SPECIFICATION
-------------------------------
*/

type dataframe [][]string

const maxCapacity = 8192 * 8192

/*
-------------------------------
MOTIF PROCESSING
-------------------------------
*/

/*
FindMotifPosOneChrom takes in a motif string and a reference string and finds all the positions for the motif
in the reference, returned as a 2D slice of ints, with each subslice being 2 long containing [start (inclusive), end (exclusive)].
In the context of this program, this means the reference is a chromosome string.
*/
func FindMotifPosOneChrom(motif string, chrom string) [][]int {
	re := regexp.MustCompile(motif)
	motif_pos := re.FindAllStringIndex(chrom, -1)

	return motif_pos
}

/*
FindArtifactsOnePos takes in a 2-length position slice and a dataframe of artifacts, as well as an integer
to designate the distance in characters to examine from the motif staring and ending positions.
It returns a sub-dataframe of all the rows that are within the motif's starting positions.
With how this function is structured, it will ONLY check the closest ends, so
we compare artifact ends to the start position, and artifact starts to the end position.
*/
func FindArtifactsOnePos(positions []int, df dataframe, within int) dataframe {
	var sub_df dataframe

	for i := range df {
		//assumes structure of the dataframe is a GTF file, meaning that start is index 3 and end is index 4
		currentRow := df[i]

		if currentRow[3] == "" {
			continue
		}

		currentStart := ConvertStringToInt(currentRow[3])
		currentEnd := ConvertStringToInt(currentRow[4])

		if positions[0] > currentEnd && (positions[0]-within) < currentEnd {
			sub_df = append(sub_df, currentRow)
		} else if positions[1] < currentStart && (positions[1]+within) > currentStart {
			sub_df = append(sub_df, currentRow)
		}
	}

	return sub_df
}

/*
FindArtifacts takes in a 2D slice of integers from FindMotifPosOneChrom, a GTF dataframe, a string to denote the chromosome,
an int radius, and a number of processors for parallel processing, returning the full dataframe of artifacts in proximity to the positions.
*/
func FindArtifactsOneChrom(positions [][]int, df dataframe, chrom string, within int, numProcs int) dataframe {
	var artifacts dataframe

	filtered_df := df.FilterGTFByChrom(chrom)

	c := make(chan dataframe)
	n := len(positions)
	for i := 0; i < numProcs; i++ {
		start := (i * n) / numProcs
		end := ((i + 1) * n) / numProcs

		go FindArtifactsParallel(positions[start:end], filtered_df, within, c)
	}

	for i := 0; i < numProcs; i++ {
		var sub_df dataframe
		sub_df = <-c
		artifacts = append(artifacts, sub_df...)
	}

	return artifacts.RemoveDuplicates() //will only remove a few, if any
}

/*
FindArtifactsParallel takes the same arguments as FindArtifacts but with slight variations. Instead of the full 2D slice
of positions it takes a subslice, and it takes a channel to put the end dataframe result.
*/
func FindArtifactsParallel(positions [][]int, df dataframe, within int, c chan dataframe) {
	var artifacts dataframe

	for i := range positions {
		artifactsOnePos := FindArtifactsOnePos(positions[i], df, within)
		artifacts = append(artifacts, artifactsOnePos...)
	}

	artifacts = artifacts.RemoveDuplicates()

	c <- artifacts
}

/*
FindGOTerms is a method for a dataframe of GTF artifacts, presumably from FindArtifactsOneChrom() or similar. It takes
a dataframe of GO terms, presumably from ReadGoTerms() and a number of processors to use, and uses parallel processing
to find the terms associated with each of the artifacts.
*/
func (df dataframe) FindGOTerms(goterms dataframe, numProcs int) dataframe {
	var results dataframe

	c := make(chan dataframe)
	n := len(df)
	for i := 0; i < numProcs; i++ {
		start := (i * n) / numProcs
		end := ((i + 1) * n) / numProcs

		go FindGOTermsParallel(df[start:end], goterms, c)
	}

	for i := 0; i < numProcs; i++ {
		var sub_df dataframe
		sub_df = <-c
		results = append(results, sub_df...)
	}

	return results.RemoveGODuplicates()
}

/*
FindGOTermsParallel is a helper fucntion for FindGOTerms() that takes in a sub-dataframe of GTF artifacts,
the whole GOterms dataframe, and a channel to put the results in afterwards.
*/
func FindGOTermsParallel(df dataframe, goterms dataframe, c chan dataframe) {
	var sub_df dataframe

	for i := range df {
		currentRow := df[i]
		currentID := currentRow[8][42:57]

		for j := range goterms {
			currentGOEntry := goterms[j]
			if currentID == currentGOEntry[0] {
				sub_df = append(sub_df, currentGOEntry)
			}
		}
	}

	c <- sub_df
}

/*
CountGOTerms is a serial method for a GOTerm dataframe (presumably from FindGOTerms) that counts the number
of times each term shows up in the dataframe. It returns a map of strings (terms) to ints.
*/
func (df dataframe) CountGOTerms() map[string]int {
	counts := make(map[string]int)

	for i := range df {
		current_term := df[i][3]
		//does the key exist in the map?
		if _, exists := counts[current_term]; !exists {
			counts[current_term] = 1
		} else {
			counts[current_term] += 1
		}
	}

	return counts
}

/*
FindTermRatios is a dataframe method for the counts of each GO term and finds their ratio out of the total number of
*/
func (df dataframe) FindTermRatios() []float64 {
	ratios := make([]float64, len(df)-1)

	for i := 1; i < len(df); i++ {
		currentRow := df[i]
		currentCount := ConvertStringToFloat(currentRow[1])
		ratios[i-1] = currentCount
	}

	totalCount := FindColTotal(ratios)
	for i := range ratios {
		ratios[i] = ratios[i] / totalCount
	}

	return ratios
}

/*
TestTermRatios is a function that takes in the counts dataframe and the background dataframe. For every term in the
counts dataframe, it finds the corresponding term in the background and runs a binomial test on the two proportions and count numbers.
It returns a slice of float slices, where the 0th slice contains z-scores and the 1st slice contains unadjusted p-values.
*/
func TestTermRatios(counts dataframe, bg dataframe) [][]float64 {
	results_full := make([][]float64, 2)

	z_scores := make([]float64, len(counts)-1)
	p_vals := make([]float64, len(counts)-1)

	raw_counts := counts.GetColumn(1)
	bg_counts := bg.GetColumn(1)

	raw_count_total := FindColTotal(MakeNumeric(raw_counts[1:]))
	bg_count_total := FindColTotal(MakeNumeric(bg_counts[1:]))

	for i := 1; i < len(counts); i++ {
		currentTerm := counts[i][0] //assumes term name is in the first column
		termIndex := bg.FindTermRowIndex(currentTerm)
		if termIndex == -1 {
			panic("Term not found in background GO terms somehow")
		}

		z := BinomZScore(counts[i], bg[termIndex], int(raw_count_total), int(bg_count_total))
		p := 2.0 * (1.0 - StandardNormCDF(z))

		z_scores[i-1] = z
		p_vals[i-1] = p
	}

	results_full[0] = z_scores
	results_full[1] = p_vals

	return results_full
}

/*
StandardNormCDF computes the standard normal cumulative distribution function.
In the context of this project, it is used for finding the p-value for a binomial test,
and the p-value result will end up being 1 minus this amount.
*/
func StandardNormCDF(x float64) float64 {
	if x > 0 {
		x = x * -1 //make x negative if it is not already
	}
	return 0.5 * math.Erfc(x/math.Sqrt2)
}

/*
BinomZScore takes in two rows from count dataframes and two integer count totals. It performs a standard
z-test on the two proportions and returns the resulting z-score.
*/
func BinomZScore(r1 []string, r2 []string, n1 int, n2 int) float64 {
	m1 := ConvertStringToInt(r1[1])
	m2 := ConvertStringToInt(r2[1])

	p1 := ConvertStringToFloat(r1[2])
	p2 := ConvertStringToFloat(r2[2])

	p_c := float64(m1+m2) / float64(n1+n2)

	num := p1 - p2
	denom := math.Sqrt((p_c) * (1.0 - p_c) * ((1.0 / float64(n1)) + (1.0 / float64(n2))))
	return num / denom
}

/*
For doing p-value correction, we're doing Holm-Bonfferoni correction, which requires us to
order our p-values to assign ranks to them. To ensure the original slice of p-values does not get
reordered, we have a function that will make a new list of the same to be operated on.
*/
func Copy(list []float64) []float64 {
	newList := make([]float64, 0)
	for _, i := range list {
		newList = append(newList, i)
	}

	return newList
}

/*
MergeSort takes a slice of floats and returns that same slice of integers, but sorted.
It's just standard mergesort.
*/
func MergeSort(inList []float64) []float64 {
	//first, get length of list
	n := len(inList)

	//check if length is 1: if length is 1, we're already sorted.
	if n == 1 {
		return inList
	}

	//otherwise, we divide this up into parts
	a := MergeSort(inList[:(n / 2)])
	b := MergeSort(inList[(n / 2):])

	return Merge(a, b)
}

/*
Merge takes two sorted lists and returns a single list sorting the elements of both lists.
Meant primarily to be used as a helper for MergeSort.
*/
func Merge(list1, list2 []float64) []float64 {
	outList := make([]float64, 0)

	//first check if both lists have elements remaining
	for len(list1) > 0 && len(list2) > 0 {
		//fetch the first element (smallest) from each list
		x, y := list1[0], list2[0]

		//compare the two elements and append the smaller one to outList
		if x < y {
			outList = append(outList, x)
			list1 = list1[1:]
		} else {
			outList = append(outList, y)
			list2 = list2[1:]
		}
	}

	//at this point, one of our lists is completely empty
	//we now append the other list to outList
	if len(list1) > 0 {
		outList = append(outList, list1...)
	}
	if len(list2) > 0 {
		outList = append(outList, list2...)
	}

	return outList
}

/*
CorrectPVals takes in a slice of n floats (p-values) and runs Holm-Bonferroni correction on the values.
It sorts them and uses the positions to apply ranks (k), then multiplies each p-value by by
(n + 1 - k) and returns this as the corrected p-value, unless it is greater than 1, in which case it returns 1.
Note that this program does NOT use Benjamini-Hochberg correction as that correction
assumes independence for each test; here, since the existence of a GO term artifact near a motif
can often exclude the existence of another artifact near that motif, this program does not assume
independent tests. Furthermore, since our P-value correction relies on finding the first instance of
a p-value in the ascending sorted list, correction tends to be more conservative.
*/
func CorrectPVals(pvals []float64, prec float64) []float64 {
	n := len(pvals)
	pvals_sorted := Copy(pvals)
	pvals_sorted = MergeSort(pvals_sorted)

	padj := make([]float64, n)
	for i := range pvals {
		p_current := pvals[i]

		rank := FindValueIndex(pvals_sorted, p_current, prec)
		if rank == -1 {
			panic("P-value not found; you may need to adjust your parameters.")
		}
		p_new := p_current * (float64(n) + 1.0 - float64(rank))
		if p_new > 1.0 {
			p_new = 1.0
		}
		padj[i] = p_new
	}

	return padj
}

/*
-------------------------------
DATA PROCESSING AND HANDLING
-------------------------------
*/

/*
While you could use packages to do this, GetKeys allows us to use base Golang to get the keys for the FASTA map.
*/
func GetKeys(m map[string]string) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}

	return keys
}

/*
ReadFasta takes in a filename for a FASTA file and returns a fasta (a map from strings to strings)
taken from the file.
*/
func ReadFasta(filename string) map[string]string {
	f := make(map[string]string)

	//standard file reading code
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer file.Close()

	index := 0
	descriptions := make([]string, 0)
	seqs := make([]string, 0)

	scanner := bufio.NewScanner(file)
	buf := make([]byte, maxCapacity)
	scanner.Buffer(buf, maxCapacity)

	for scanner.Scan() {
		currentRow := scanner.Text()
		if index%2 == 0 {
			descriptions = append(descriptions, currentRow)
		} else {
			seqs = append(seqs, currentRow)
		}
		index++
	}

	if len(descriptions) != len(seqs) {
		panic("Inputted FASTA file does not have equal number of descriptions as sequences.")
	}

	for i := 0; i < len(descriptions); i++ {
		f[descriptions[i]] = seqs[i]
	}

	return f
}

/*
ReadGTF takes in a filename for a string and returns a 2-dimenstional array of strings to store the information
from the file. With how it is formatted, it will automatically create heading labels corresponding to a standard GTF
(["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]).
*/
func ReadGTF(filename string) dataframe {
	var df dataframe
	currentIndex := 0

	//standard file reading code
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer file.Close()

	colNames := make([]string, 9)
	colNames[0] = "chr"
	colNames[1] = "source"
	colNames[2] = "feature"
	colNames[3] = "start"
	colNames[4] = "end"
	colNames[5] = "score"
	colNames[6] = "strand"
	colNames[7] = "frame"
	colNames[8] = "attribute"
	df = append(df, colNames)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		currentRow := scanner.Text()
		entries := strings.Split(currentRow, "\t")

		//make a new row
		newRow := make([]string, len(entries))
		df = append(df, newRow)
		for i := range entries {
			df[currentIndex][i] = entries[i]
		}

		currentIndex++
	}

	err2 := scanner.Err()
	if err2 != nil {
		panic(err2)
	}

	return df
}

/*
ReadGOTerms takes in a filename for a string and returns a 2-dimenstional array of strings to store the information
from the file. With how it is formatted, it will automatically create heading labels corresponding to a standard GO
term TSV (["id", "type", "term_id", "term_name", "description", "source"]).
*/
func ReadGOTerms(filename string) dataframe {
	var df dataframe
	currentIndex := 0

	//standard file reading code
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer file.Close()

	colNames := make([]string, 6)
	colNames[0] = "id"
	colNames[1] = "type"
	colNames[2] = "term_id"
	colNames[3] = "term_name"
	colNames[4] = "description"
	colNames[5] = "source"
	df = append(df, colNames)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		currentRow := scanner.Text()
		entries := strings.Split(currentRow, "\t")

		//make a new row
		newRow := make([]string, len(entries))
		df = append(df, newRow)
		for i := range entries {
			df[currentIndex][i] = entries[i]
		}

		currentIndex++
	}

	err2 := scanner.Err()
	if err2 != nil {
		panic(err2)
	}

	return df
}

/*
ReadCSV takes in a string name that will be the filename for the CSV we want to read.
It returns a dataframe where each sub-array in the array corresponds to a row of the CSV.
*/
func ReadCSV(filename string) dataframe {
	var df dataframe
	currentIndex := 0

	//standard file reading code
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		currentRow := scanner.Text()
		entries := strings.Split(currentRow, ",")

		//make a new row
		newRow := make([]string, len(entries))
		df = append(df, newRow)
		for i := range entries {
			df[currentIndex][i] = entries[i]
		}

		currentIndex++
	}

	err2 := scanner.Err()
	if err2 != nil {
		panic(err2)
	}

	return df
}

/*
ReadMotifsFromFile is a function that takes in a file of motifs, one motif per line.
It returns the capitalized motifs in a slice of strings. It also strips the lines and checks if the length is 0.
*/
func ReadMotifsFromFile(filename string) []string {
	motifs := make([]string, 0)
	//standard file reading code
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		currentLine := scanner.Text()
		currentMotif := strings.TrimSpace(currentLine)
		currentMotif = strings.ToUpper(currentMotif)
		motifs = append(motifs, currentMotif)
	}
	return motifs
}

/*
RemoveDuplicates is a dataframe method that removes rows that are completely duplicated.
For removing rows based on specific column index values, use a different function.
*/
func (df dataframe) RemoveDuplicates() dataframe {
	var sub_df dataframe

	for i := range df {
		currentRow := df[i]
		if i == 0 {
			sub_df = append(sub_df, currentRow)
			continue
		}

		hasEqual := false
		for j := range sub_df {
			if AreSlicesEqual(currentRow, sub_df[j]) {
				hasEqual = true
				break
			}
		}
		if !hasEqual {
			sub_df = append(sub_df, currentRow)
		}
	}

	return sub_df
}

/*
The RemoveGODuplicates() method is designed for a GO term dataframe, in that the duplicate checking is done only for columns
0 (GeneID) and 2 (GOTermID). Typically when these are the same, the only difference (if any) is in the source (InterPro/OMA.)
*/
func (df dataframe) RemoveGODuplicates() dataframe {
	var sub_df dataframe

	for i := range df {
		currentRow := df[i]
		if i == 0 {
			sub_df = append(sub_df, currentRow)
			continue
		}

		hasEqual := false
		for j := range sub_df {
			currentCheckRow := sub_df[j]
			if currentRow[0] == currentCheckRow[0] && currentRow[2] == currentCheckRow[2] {
				hasEqual = true
			}
		}

		if !hasEqual {
			sub_df = append(sub_df, currentRow)
		}
	}

	return sub_df
}

/*
Helper function for RemoveDuplicates() that checks if any two slices are equal in their entirety.
*/
func AreSlicesEqual(s1 []string, s2 []string) bool {
	if len(s1) != len(s2) {
		return false
	}
	for i := range s1 {
		if s1[i] != s2[i] {
			return false
		}
	}
	return true
}

/*
FilterGTFByChrom is a dataframe method that takes in a chromosome string.
It returns a dataframe with only rows containing that string in the first column.
*/
func (df dataframe) FilterGTFByChrom(chrom string) dataframe {
	var sub_df dataframe

	for i := range df {
		currentRow := df[i]
		if currentRow[0] == chrom {
			sub_df = append(sub_df, currentRow)
		}
	}

	return sub_df
}

/*
WriteCSV takes in a string dataframe and a file name and writes a CSV file of that name. It's really that simple.
*/
func WriteCSV(df dataframe, filename string) {
	file, err := os.Create(filename)
	if err != nil {
		panic("Error in creating CSV of specified filename.")
	}

	writer := bufio.NewWriter(file)

	for i := range df {
		for j := range df[i] {
			n := len(df[i]) - 1
			fmt.Fprint(writer, df[i][j])
			if j != n {
				fmt.Fprint(writer, ",")
			}
		}
		fmt.Fprintln(writer, "")
	}

	writer.Flush()

	file.Close()
}

/*
WriteCounts() is a function that takes a map of counts (presumably from GetCounts) and a filename, makes a dataframe out
of the map, and sends that dataframe to be written with that filename as a CSV.
*/
func WriteCounts(counts map[string]int, filename string) {
	df := CountsToDF(counts)
	WriteCSV(df, filename)
}

/*
CountsToDF is a function that takes a map of counts (presumably from GetCounts and a filename, makes a
dataframe out of it with the first column being the term name and the second the count, and returns
the counts as a dataframe.
*/
func CountsToDF(counts map[string]int) dataframe {
	var df dataframe

	columnNames := make([]string, 2)
	columnNames[0] = "name"
	columnNames[1] = "count"
	df = append(df, columnNames)

	//since our GetKeys function takes maps of strings to strings, it doesn't work for our counts mapping
	keys := make([]string, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}

	for k := range keys {
		newRow := make([]string, 2)

		//terms might have commas, which we need to avoid
		currentTerm := keys[k]
		containsComma := strings.ContainsRune(currentTerm, ',')
		if containsComma {
			currentTerm = strings.ReplaceAll(currentTerm, ",", ";")
		}
		newRow[0] = currentTerm

		c := strconv.Itoa(counts[keys[k]])
		newRow[1] = c

		df = append(df, newRow)
	}

	return df
}

/*
GetColumn is a method for a dataframe that takes in an integer index and returns a column corresponding to that
column position in the dataframe.
*/
func (df dataframe) GetColumn(index int) []string {
	//error statements.
	if len(df) < 1 {
		panic("Empty dataframe passed to GetColumn().")
	} else if len(df[0]) < 1 {
		panic("Empty dataframe passed to GetColumn().")
	} else if len(df[0]) < index {
		panic("Index too big for GetColumn().")
	}

	// ok now we can do stuff
	col := make([]string, len(df))
	for i := range df {
		col[i] = df[i][index]
	}

	return col
}

/*
NumCols is a method for a dataframe that returns the number of columns.
*/
func (df dataframe) NumCols() int {
	if len(df) < 1 {
		panic("error: empty dataframe given to NumCols()")
	}
	return len(df[0])
}

/*
NumRows is a method for a dataframe that returns the number of rows.
*/
func (df dataframe) NumRows() int {
	return len(df)
}

/*
RemoveNAs is a dataframe method that takes in an integer index and removes all rows where we have an
NA in that index. For the context of GTFs and GO terms files, this is unnecessary as these files
are typically fully populated. NAs can be interpreted as the actual string "NA" or an empty string.
*/
func (df dataframe) RemoveNAs(index int) dataframe {
	//error statements.
	if len(df) < 1 {
		panic("Empty dataframe passed to RemoveNAs().")
	} else if len(df[0]) < 1 {
		panic("Empty dataframe passed to RemoveNAs().")
	} else if len(df[0]) < index {
		panic("Index too big for RemoveNAs().")
	}

	//ok now we can actually do stuff
	var newdf dataframe
	//range through and add valid rows
	for i := range df {
		if df[i][index] != "" && df[i][index] != "NA" {
			newdf = append(newdf, df[i])
		}
	}

	return newdf
}

/*
DropColumn is a function that takes in a dataframe and an index for a column of the frame.
It returns all columns but that one.
*/
func (df dataframe) DropColumn(index int) dataframe {
	var newdf dataframe

	//range through rows
	for i := range df {
		var newRow []string

		//range through cols
		for j := range df[i] {
			if j != index {
				newRow = append(newRow, df[i][j])
			}
		}

		newdf = append(newdf, newRow)
	}

	return newdf
}

/*
Select is a function that is basically the opposite of DropColumn; it takes an array of names and selects only those columns.
*/
func (df dataframe) Select(names []string) dataframe {
	var newCols dataframe

	//range through and get all columns
	for i := range names {
		if ColIn(names[i], df) {
			newCol := df.GetColumn(df.GetIndex(names[i]))
			newCols = append(newCols, newCol)
		}
	}

	newDF := JoinCols(newCols)
	return newDF
}

/*
ColIn is a function that returns true if an input string is in the input dataframe's first row.
Otherwise, it returns false.
*/
func ColIn(name string, df dataframe) bool {
	//range through all col names and compare to inpuit
	for i := range df[0] {
		if df[0][i] == name {
			return true
		}
	}

	//else, we haven't found the name, so return false
	return false
}

/*
JoinCols is a function that takes the columns given by Select() and GetColumn() and joins them
into an appropriate dataframe.
*/
func JoinCols(cols1 dataframe) dataframe {
	var newDF dataframe

	//range through and append as necessary
	for i := range cols1[0] { //first column given, which indexes the rows. Does that make sense? Kindof?
		var newRow []string

		for j := range cols1 { //range through the columns now
			newRow = append(newRow, cols1[i][j])
		}

		newDF = append(newDF, newRow)
	}

	return newDF
}

/*
GetIndex is a function that takes in a string value of a name for a column in a dataframe and
returns the integer index of that variable.
*/
func (df dataframe) GetIndex(name string) int {
	for i := range df[0] {
		if df[0][i] == name {
			return i
		}
	}

	panic("Error: name does not exist in dataframe given to GetIndex().")
	return -1
}

/*
IsNumeric is a function that takes in a slice of strings and returns whether or not
one of them is a number or not, by using the error of ParseFloat.
*/
func IsNumeric(col []string) bool {
	for i := range col {
		_, err := strconv.ParseFloat(col[i], 64)
		if err != nil {
			return false
		}
	}

	return true
}

/*
MakeNumeric is a function that takes in a slice of strings and returns a slice of those
strings converted to floats.
*/
func MakeNumeric(col []string) []float64 {
	newCol := make([]float64, len(col))

	for i := range col {
		newVal := ConvertStringToFloat(col[i])
		newCol[i] = newVal
	}

	return newCol
}

/*
FindColTotal is a function that takes a slice of floats and returns the total.
*/
func FindColTotal(nums []float64) float64 {
	total := 0.0
	for i := range nums {
		total += nums[i]
	}
	return total
}

/*
AddNumericColumn is a dataframe method that takes a slice of floats and a string column name. It returns the dataframe
with the numeric column converted to strings and appended at the end of the dataframe.
*/
func (df dataframe) AddNumericColumn(newCol []float64, name string) dataframe {
	var df_new dataframe //copy over everything to prevent overlapping, in case
	df_new = append(df_new, df...)

	df_new[0] = append(df_new[0], name)
	for i := 1; i < len(df_new); i++ {
		currentVal := strconv.FormatFloat(newCol[i-1], 'f', -1, 64)
		df_new[i] = append(df_new[i], currentVal)
	}

	return df_new
}

/*
FindTermRowIndex is a dataframe method that takes in a for a GO term and returns the row index of the
term in the dataframe. It returns -1 if the term is not found.
*/
func (df dataframe) FindTermRowIndex(term string) int {
	for i := 1; i < len(df); i++ {
		//assumes the term is in the first column
		if term == df[i][0] {
			return i
		}
	}
	return -1
}

/*
FindValueIndex is a function that takes in a slice of floats and a singular float
for which you want to find its position. It returns that position, returning -1 if it is not found.
Note that for this, the user specifies a precision, which for now is 0.00001.
*/
func FindValueIndex(values []float64, term float64, prec float64) int {
	for i := 1; i < len(values); i++ {
		//assumes the term is in the first column
		diff := math.Abs(term - values[i])
		if diff < prec {
			return i
		}
	}
	return -1
}

/*
----------------------------
MISCELLANEOUS FILE AND STRING HANDLING
----------------------------
*/

/*
DoesFileExist() checks to see if a file exists. It panics if there is an error checking the file, such as
with permissions being denied.
*/
func DoesFileExist(filePath string) bool {
	_, err := os.Stat(filePath)
	if err == nil {
		return true
	} else if errors.Is(err, os.ErrNotExist) {
		return false
	} else {
		// Another error occurred (e.g., permission denied)
		panic("Error checking file: %v")
	}
}

/*
ConvertStringToInt is a simple function that takes in a string and returns the integer it converts to.
It's essentially a wrapper for strconv.Atoi(s), with appropriate error handling.
*/
func ConvertStringToInt(s string) int {
	i, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return i
}

/*
ConvertStringToFloat is a simple function that takes in a string and returns the float it converts to.
It's essentially a wrapper for strconv.ParseFloat(), with error handling.
*/
func ConvertStringToFloat(s string) float64 {
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic(err)
	}
	return f
}

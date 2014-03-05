Centromere Repeat Evolution in the Grasses

October 2013-February 2014
by: Paul Bilinski

Below is a running log/journal of the final analyses performed for the Hufford et al
centromere evolution paper.  Entries are not made in chronological order, but rather
broken up by tasks completed.  Different form of organization than the Waters et al
paper, but hopefully better.

Searching the document:

--- indicates new task
!!! indicates pressing question to be addressed, probably by jeff
### indicates broad project segment

###Getting Genomic Abundance

--- Abundance in the different Centromeric Retros

Initial examination of the CRM in rice and sorghum showed that there is not much difference
in using either rice or maize or sorghum (see CRR CRS CRM testing in PB farm dir).  To
test this on empirical data, I ran each of our files (first iteration, obviously flawed
as a lot of the lanes did not work) against the different homologous regions of the TEs.

###Building De novo Sequence

--- De Novo Building with Mira

Submit script via sbatch Submit_assembly.sh

	#!/bin/bash
	#SBATCH -D /home/pbilinsk/huffwork/Assembly
	#SBATCH -o /home/pbilinsk/huffwork/Assembly/slurmlog/assembly-stout-%j.txt 
	#SBATCH -e /home/pbilinsk/huffwork/Assembly/slurmlog/assembly-sterr-%j.txt
	#SBATCH -J Assemble
	#SBATCH -p serial
	set -e
	set -u


	mira manifest.conf >&log_assembly.txt
	
The manifest file looks like:

	project = RimmarepeatASS 
	job = genome,denovo,accurate
	parameters = -highlyrepetitive -NW:cnfs=no -NW:mrnl=200 -HS:mnr=no

	readgroup = Rimma
	data = /home/pbilinsk/huffwork/Data/FWD/MH_RIMMA0019.J.fastq
	technology = solexa	

This produces a folder that contains project_results.

This was repeated for each of the assemblies, with suitable paths.  All files are on farm,
in the directory huffwork under user pbilinsk

--- Using Tandem repeat finder on the assembly file from mira and filter

In the _results folder, you will find the outputs of the assembly.  One of them will be
the unpadded file.  example:

RimmarepeatASS_out.unpadded.fasta
ZperRepeat_out.unpadded.fasta
ApludaRepeat_out.unpadded.fasta
SorghumRepeat_out.unpadded.fasta
TperuRepeat_out.unpadded.fasta
PhyloRepeat_out.unpadded.fasta
TdactRepeat_out.unpadded.fasta
TriturRepeat_out.unpadded.fasta

Download this file from the cluster, and run trf on it locally.

./trf407b.macos64 RimmarepeatASS_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 ZperRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 ApludaRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 SorghumRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 TanderRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 TperuRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 PhyloRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 TdactRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 TriturRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h


Next, take the .dat produced by TRF and filter it based on the length requirements we see
in Melters et al 2013.  They found that the shortest tandem repeat in a plant was 40bp
long, in ricin.  The script TRF_parser.pl will extract all sequences whose tandem repeat
length is more than 40bp.  

perl TRF_parser.pl RimmarepeatASS_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Rimmarepeat_TRFfinds.fasta
perl TRF_parser.pl ZperRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Zperrepeat_TRFfinds.fasta
perl TRF_parser.pl ApludaRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Apludarepeat_TRFfinds.fasta
perl TRF_parser.pl SorghumRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Sorghumrepeat_TRFfinds.fasta
perl TRF_parser.pl TperuRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Tperurepeat_TRFfinds.fasta
perl TRF_parser.pl TanderRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Tanderrepeat_TRFfinds.fasta
perl TRF_parser.pl TdactRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Tdactrepeat_TRFfinds.fasta
perl TRF_parser.pl TriturRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Triturrepeat_TRFfinds.fasta

perl TRF_parser.pl PhyloRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Phylorepeat_TRFfinds.fasta


This will produce a file where each line is a tandem repeat and its name is >textassembly#

--- Abundance of non-Knob elements, what is in the top 10?

First, we want to understand what the most abundant of the tandem repeats is.  By far, the
most common element with be the knob repeat, and therefore we have to mask out this sequence.  
To do this, we use BLAST and the 34 annotated knob sequences from NCBI.  

makeblastdb -in 34KnobsonNCBI_renamed.txt -dbtype 'nucl' -parse_seqids

With the DB made, I then blast our tandem repeats against the database:

blastn -query Rimmarepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_rimmatrf
blastn -query Zperrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_zpertrf
blastn -query Apludarepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_apludatrf
blastn -query Sorghumrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_sorghumtrf
blastn -query Tanderrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_tandertrf
blastn -query Tperurepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_tperutrf
blastn -query Tdactrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_tdacttrf
blastn -query Triturrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_triturtrf

blastn -query Phylorepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_phylotrf

Transfer files to my computer, the playdir directory, and run: 

perl Blast_DBparser.pl DB_knobs_vs_rimmatrf > assemblynamesinDB_knobs_vs_rimmatrf.txt
perl Blast_DBparser.pl DB_knobs_vs_zpertrf > assemblynamesinDB_knobs_vs_zpertrf.txt
perl Blast_DBparser.pl DB_knobs_vs_apludatrf > assemblynamesinDB_knobs_vs_apludatrf.txt
perl Blast_DBparser.pl DB_knobs_vs_sorghumtrf > assemblynamesinDB_knobs_vs_sorghumtrf.txt
perl Blast_DBparser.pl DB_knobs_vs_tandertrf > assemblynamesinDB_knobs_vs_tandertrf.txt
perl Blast_DBparser.pl DB_knobs_vs_tperutrf > assemblynamesinDB_knobs_vs_tperutrf.txt
perl Blast_DBparser.pl DB_knobs_vs_tdacttrf > assemblynamesinDB_knobs_vs_tdacttrf.txt
perl Blast_DBparser.pl DB_knobs_vs_triturtrf > assemblynamesinDB_knobs_vs_triturtrf.txt

perl Blast_DBparser.pl DB_knobs_vs_phylotrf > assemblynamesinDB_knobs_vs_phylotrf.txt

Notes: tdact has no hits.  sorghum has 1 hit.  The rest have a ton.  Those 2 were done by
hand for ReadyForMosaik files.  Second note, phylo has no knobs.

I print out the names of each sequence that has a hit longer than 30bp to one of the 
repeats.  This is printed out, duplicates are removed via textwrangler (process duplicate
lines, delete duplicate lines, print to new file named assemblynamesinDB_GRAMvsRIMMA_nodupl.txt).  ###DID NOT RENAME
The original names of all of the assemblies are taken using 

grep ">" Rimmarepeat_TRFfinds.fasta > Rimma_assemblynames_forsubset.txt
grep ">" Zperrepeat_TRFfinds.fasta > Zper_assemblynames_forsubset.txt
grep ">" Apludarepeat_TRFfinds.fasta > Apluda_assemblynames_forsubset.txt
XXgrep ">" Sorghumrepeat_TRFfinds.fasta > Sorghum_assemblynames_forsubset.txtXX
grep ">" Tanderrepeat_TRFfinds.fasta > Tander_assemblynames_forsubset.txt
grep ">" Tperurepeat_TRFfinds.fasta > Tperu_assemblynames_forsubset.txt
XXgrep ">" Tdactrepeat_TRFfinds.fasta > Tdact_assemblynames_forsubset.txtXX
grep ">" Triturrepeat_TRFfinds.fasta > Tritur_assemblynames_forsubset.txt

Don't forget to get rid of the > in text wrangler.

The list is then imported into R, where I load 
the list of names of all of the test assemblies.  Use the short script:

	Orig <- read.csv("Rimma_assemblynames_forsubset.txt", header=FALSE)
	remain <- read.csv("assemblynamesinDB_knobs_vs_rimma_nodupl.txt", header=FALSE)
	notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
	write.csv(notknownrepeats,file = "nonrepeatsknob.csv")

This will create a csv file with additional columns and such.  I just created a new file
from these names by copy pasting all of them into NAME_nonknobassemblynames.txt, and adding the > 
with replace of testassemblies for >testassemblies.  Now, we need to use those names that
dont have a hit against all the repeats to get the assemblies are unknown.  We have the
script Fetch_nonrepeatassemblies.pl.

perl Fetch_nonrepeatassemblies.pl Rimma_nonknobassemblynames.txt Rimmarepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl Zper_nonknobassemblynames.txt Zperrepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl Apluda_nonknobassemblynames.txt Apludarepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl Sorghum_nonknobassemblynames.txt Sorghumrepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl Tander_nonknobassemblynames.txt Tanderrepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl Tperu_nonknobassemblynames.txt Tperurepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl Tdact_nonknobassemblynames.txt Tdactrepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl Tritur_nonknobassemblynames.txt Triturrepeat_TRFfinds.fasta > test.txt

Then, test whether those names are the same in each of the files.

grep ">" test.txt > compareafter
grep ">" Tritur_nonknobassemblynames.txt > comparebefore
diff comparebefore compareafter 

In test.txt, I found that the very last entry was wrong, and replaced it... Will have to do this for
each of the individuals.  I replaced:

>testassemblies0
CTGGCTCGGGCCGATTCCAGCGTAAACCGTGAGCTAAAACAGCGTAAATGAGTATAGAAATTTAGCGTAAATCTTATATCTGTTTTGTAACAGCAAATGAGGCCTAAAATTACGGCGTGAAAATTGTGTGCAGCCGATCGTGCACGGGTC

>testassemblies7419
CTTCGCCTTATATTTCGGCGGAATCAGCGTTGACTTTTCGCG


mv test.txt ReadyForMosaik_rimma.txt
mv test.txt ReadyForMosaik_zper.txt
mv test.txt ReadyForMosaik_apluda.txt
mv test.txt ReadyForMosaik_tander.txt
mv test.txt ReadyForMosaik_tperu.txt
mv test.txt ReadyForMosaik_tritur.txt

#Sorghum was hand fixed because only 1 sequence needed removal, and tdact had none. Makes 8.
#Phylo had none, so that makes 9.

Renamed the file since it is now ready for mosaik!  Moved onto the cluster, and mapped.

###Abundance of the De novo sequence 

---  Identifying the most common tandem repeat via mosaik

After renaming, moved to the cluster and made it into mosaik ready .dat.  Operations were
executed in the ~/huffwork/Mosaikmapping/ directory on the cluster.  Then build the 
sequences and submit the script via (all files must go to correct dir):

	./MosaikBuild -fr ReadyForMosaik_NAME.txt -oa ReadyForMosaik_NAME.dat
	sbatch Submitfulltest.sh ReadyForMosaik_nonknob.dat

The references are in the ~/huffwork/References directory.  Once mapped, in each separate folder the
alignments produce a *.loc file.  This will give a listing of each assembly and how many reads
mapped to it.  Do a head ___ > RanksNAME.csv  for each of the loc files.  Download them to
a local directory, and open them up in text wrangler.  Format this so you can open it as a
csv in excel.  Order them, find the top hit, and locate it in the TRF output file.  Find
the monomer version of the longer hit, those are shown below.

Apluda: like sorghum
>testassemblies383(long)
GAAACTCATTTCGGCCTGTTTGGAGACTCAGATTTATCCCAGTGCAACATAGGTGCACGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGTGGAAACTAGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCAGTGCAACATAGGTGCCCGGTTTGCATCGAATGTACCATCAGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGTGGAAACTCGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCAGTGCAACATAGGTGCCCGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGTGGAAACTCGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCGGTGCAACATAGGTGCCCGATTTGCGTCCAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTACTAACGCGGGTCCAGTGGAAACTCGTTTTGGCCTGTTTGGAGACTCAGATTTATCCCGGTGCAACATAGGTGCCCGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGT
>testassemblies384(short, 138bp)
TGGAAACTCGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCGGTGCAACATAGGTGCCCGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAG

Phylo: blasts to bamboo repetitive region, not in melters, not annotated as anything
>testassemblies114(long)
TGTGAAACTTTCCAGAACGACGCGTACCCTAAGCTCCAATTTATGAAATCGAGGCACAGAAAGAAAATCACACCAAATAAGTATAATGGACTTCTGGTGTGCCCTACGTCATTTCTAGGTGCATTTGGAATAATTCCGAGCAAGTGTGAAACATTCCAGAACGACGGGTACCCAGAGCTCCACTTTATGAAATCGAGGCATGGAAAGAAAATGAACATCAAACTCCCATAATGGACTTCTAGTGTGCCCTACCTCATTTCTAGGTGCATTTGGAATCTTTCCGAGCCGGTGTGAAACTTTCCAGAACGACGCGTACCCTGAGCTCCAATTTATGAAATCGAGGCACGGAAAGAAAATCAAAACCAAACACCCATAATGGACTTCTAGTGTGCCCTACATCATTTCTAGGTGCATTTGGAATAATTCCGACCAGGTGTGAAACTTTCCAGAACGACGGGTACCCAGAGCTCCACTTTTGAAATCGAGGCACGGAAAGAAAATCAACACCAAATAACCATAATGGACTTCTAGTGTGCCCTACCTCATTTCTAGGTGCATTTGAATCATTCCCAGCAAGTGTGAAACATTCCAGAACGACAGGTACCCACAGCTCCAATGTACGAAATCGAGGCACGGAAAGAAAATCAAAACCAAACACCCATAATGGACTTCTACTGTGCCCTACCTCATTTCTAGGTGCATTTGGAATAATTCCGAGCAAGTGTGAAACATTCCAGAACGACGCGTACCCAGAGCTCCAATGTATGAAATCGAGGCACGGAAAGAAAATCAAAACCAAACACCCATAATGGACTTCTAGTGTGCTCTACCTCATTTCTAGGTGCATTTGGAATAATTCTGAGCAGGTGTAAAACTTTGCAAAGCCACGGGAACCCAGAGCTCCACATTATGAAATCGAGGCACGGAAAGAAAATCAATGCAAACACCCAAAATGGACTTCTGGTGTGCCCTATGTCATTTCTAGGTGCATTTCGAATAATTCCGAGCAGGTGTGAAACTTTGCTGAACGACGGTATCCAGAGCTCCACTTTACGAAATCAAGGCACAGAAAGACAATCAATACGAAACACCCAAAATAGACTTCTAGTGTGCCCTACCTCATTTCTAGGTGCATTTGGAATCCTTCCGAGCCGG
>testassemblies112(short, 145bp)
TGTGAAACTTTCCAGAACGACGCGTACCTAAGCTCCAATTTATGAAATCGAGGCACGGAAAGAAAATCAACACCAAATACCCATAATGGACTTCTAGTGTGCCCTACGTCATTTCTAGGTGCATTTGGAATAATTCCGAGCCGG

Rimma: top is ribosomal, second is centc
>testassemblies82(Top abundance... ribosomal)
CGACGGGCGGCAGGGCCGTTCGGCCTTGCGTCTGGGCTGAAAACGAGGGGTTCGCCATGGCGCACGGGCCGAAAACGGAGGCTCGGGCACGAACTCGAAAATAAGCTAAGTCCAAGCGTGTGGAAAGACACCGAACCTAAAGTGCATGTCTTGAGTGAAGCGCAAGGTGGAACGGAGGGAAAACGGGAGGAGGCGCGCGGCGACGAGCGGCAGGGCCGTTCGGCCTTGCGCCTGGGCTGAAAACAAGGGGTTCGCCATGGCGCACGGGCCGAAAACCGAGGTTCGGGCATGGCCCCGGAAATAAGCTAAGTCCAAGCGTGTGGAAAGACACCGAAGGTAATATGTATGTCTTGGGTGAAGGGCATGGCGGAACGGAGGGAAAACGGCACGGAGGCGCAAGGCGACGGGCGGCATGGCTGTTCGGCCTTGCGTCTGGGCCGAAAATGAGGGGTTCGCCCATGGCGCACGGGCCGAAAACTGAGGCCCGGGCACGACCCCGAAAATAAGCTAAGTCCAAGCGTGTGGAATGGAAAGACAACGAAGCTAAGGTGTATGTAGGCTTGAGTGAAGGGCAAGGTGGAACGGAGGGAAAACGGGAGGAGGCGCGCGGCGACGGGCGGCAGGGCTGTTCGGCCTTGCGTCTGGGCTGAAAACGAGGTGTTCGCCATGGCGCACGGGCCGAAAACGGAGGCTCGGGCACGAACTCGAAAATAAGCTAAGTCCAAGCATGTGGAAAGACACCGAACATAAAGTGTATGTCTTTGGTGAAGGGCACAGTGGAACGGAGGGAAAACGACACGGAGGCACGCGACGAACGGAGGCTCGGGCACGGAGGCGCGCGGCGACGGGCGGCATGGCTGTTCGGCCTTGCGTTTGGGCTGAAAACGAGGCGTTCGCCATGGCGCGCGGGCCGAAAACAACGGTTCGGCCACGGCCTCGGAAATAAGCTAAGTCCAAGCGTGTGGAAAGACACCGAAGCTAATGTGTAAGTCTTGGGTGAAGGGCACGGTGGAACGGAGGGAAAACGACACGGAGGCACGCGA
>testassemblies302(Second, long repeat centc)
TTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGATACCCGTTTTCGGAATGGGTGACGTGCGGCAACGAAATTGCGCGAAAACAACCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACAAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCGATACCCGTTTTCGGAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTACGAAAAACGAAGAAATGGTTTCGGTGGCGAAAACTCGTGCTTTGTATGCACCCCGACACCCGTTTTCGGAATGGGTGACGTGCGACAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACAACCGTTTTCAGAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAATCCGTGCTTTGTATGCACCCGACACCCGTTTTCGGAATGGGTGACGTGCGACAACGAAATAGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTACGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCATTGTATGCACCCGACAACAGTTTTCGGAATGGGTGACGTGCGACAACGAAATTGCGAGAAACCACCACAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCAAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGAGCTTTGTATGCACCCAATACCCGTTTTCGGAATGGGTGACGTGCCGCAACGAAATTGCGCGAAAACACCACAAACATGAGTTTAGGACCTAAAGTAGTGGAGTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCATGCTTTGTATGCACCCCGATACCCGTTTTCGGAAGGGGTGACGTGCGGCAACGAAATTGCACGAAACCAACCCAAACATGAGTTTTGGACCTAAAATAGTGGAGTGGGCATGTTCG
>testassemblies301(Short, 155bp)
TTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACACCCGTTTTCGGAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCG

Sorghum: sorghum cent
>testassemblies40(long)
GGGTGCATCCAAATTGATTTCTAAGCATATGGTACGTTCCATGCAAACCGTGCACCTATCTTGCATCAAGATTAGCACTATCTCCAAACAAAGCAAACTGAGCTTCCACTTGAGCCCCTTTACCCAGGAGTATCATCGGGTGCATCCAAAATGGTTCCTTAGCCTATGATGCATTAGGCGCAAACTGTGTACCTATCTTGCACCAAAACTAACTCTGTCTCCAAACAAACCAAAGCGAGATTCCATATGACACATGTCAACTAGGAGTTCTATTGCGATGCACCCAAATAGATTTCTTAGCATATGGTATGTTCCATGAAAACCGTGCACCTATCTTGCATCAAGATTAACACTATCTCCAAAAAGAACAAACCGAGCTTCCACTTGAGCCTCTTCACCGAAGTACCAACAGGTGCGTCCAAAATGGTTTCTTAGACTATGCTGCATTAGGCGCAAACCATGCACATATCTTGCACCGAAACTAACACAATCTCCAAACAGACCAAAGCAAGATTCCATATGACACACGTCATCAAGGAGTTCCATCAGGTGCATCCAAATTGATTTCCAAGCATATGGTATGTTCCGTGCAAACCGTGCACCTATCTTGCATCAAGATTAGCACTATCTCCAAACAGACCGAACAGACCATCCACTTGAGCCCCTTCACCTAAGAATACCATCAAGTGCGTCCAAAATGGTTTCTGAGCATATGGTGCATTAGGAGCAAACGATGCACCAATCTTCCACCGAAACTAACAATGTCTCCAAACAGACCGAAGCGAGATTCCAAATGACACACATGATCAAGGAGTTCCATC
>testassemblies39(short, 137bp)
GCAAACCTGTGCACCTATCTTGCACCGAAACTAACACTGTCTCCAAACAGACCAAAGCAGATTCTATATGACACACGTCATCTAGGAGTTCCATCAAGTGCGTCCAAATTGATTTCTAAGCATATGGTATGTTCCAT

Tander: maize centc
>testassemblies157(long)
AAAAAGAAGAAATGGTTCCGGTGGCAAAAACTTATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAAGAAGAAATGGTTCCGATGGCAAAAACTCATGCCATATATGCACCCCGATACCCGTTTTCGGAATGGGTGACGTGCAGCAACGAAATGGCACGAAACCACCCCAAACATGAGTTTTGGACCCAAAGTAGTGCATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACAAAAAGCACGAAACCACCCCAAACATGAGTTTTGGACCTAAACTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGATTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCAGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCGCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCATTGCGA
>testassemblies156(short, 156bp)
AAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGA

Tdact:
>testassemblies122(343bp, only tandem repeat)
CGCTTTGAATGCTGCATACTGAACACAAAAGAAGTCCGGAGTTCAAATAAGTTTAAAAAACATTGAAGTGCCCGTGTAACAGATGAGTTCTCGTCCGAAACCCTGATACTCCGAAAGAGATTGTCCAGTTTGTACACGAAGTGCGTCCAGTTTTTGCCGTGACCCTCTCTACTCTTTCGCACATGCTATGCGGGTGAAATGATGATACCATGCCAAGTTTCAACATTTTCAGAATTCATTTTGTAGTGATTTTCAATTTCACGGTCATTTAGCTCTCTAAACAATTAGGTAAATGACCGAAAAACAACAAATGATGTCAGAACATGTTGGAAATTGATGACGT

Tperu: maize like centc
>testassemblies111(long)
AACCATTTCTTCCTTTTTCGCAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTGAAACGGGTATCGGGGTGCATACAAGGCATGAGTTTTTGCCACCGAAACCATTTCTTCCTTTTTCACAACGAACATGCCCAATCCACTACTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTCCGAAAACGGGTATCGCGGTGCATACAAGGCATGAGTTTTTGCCACCGGAACCATTTCTTCTTTTTTCACAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGCGGAGGTTTCGCACCATTTCGTTGCCGCACGTCACCCATTGCGAAAACGGGTATCGGGGTGCATACAAGGCACGAGTTTTTGCCACCGC
>testassemblies110(short, 156bp)
GGAACCATTTCTTCCTTTTTCGCAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTGCGAAAACGGGTATCGGGGTGCATACAAGGCATGAGTTTTTGCCACC

Tritur: maize centc
>testassemblies1141(long)
AACCATTTCTTCTTTTTTCGCAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTGTGAAAACGGGTATCGGGGTGCATACAAGGCATGAGTTTTTGCCACCGGAACCATTTCTTCCTTTTTCGCAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTGCGAAAACGGGTATCGGGGTGCATACAAGGCATGAGTTTTTGCCACCGGAACCATTTCTTCTTTTTCGCAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTGCGAAAACGGGTATCGGGGTGCATACAAGGCATGAGTTTTTGCCACCGGAACCATTTCTTCTTTTTTCGCAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTGCGAAAACGGGTATCGGGGTGCATACAAGGCATGAGTTTTTGCCACCGG
>testassemblies1139(short, 156bp)
AACCATTTCTTCTTTTTTCGCAACGAACATGCCCAATCCACTACTTTAGGTCCAAAACTCATGTTTGGGGTGGTTTCGCGCCATTTCGTTGCCGCACGTCACCCATTGCGAAAACGGGTATCGGGGTGCATACAAGGCATGAGTTTTTGCCACCGG

(putativetrit: CTGGCCTTGAGAAGACGTTCGAAACAAAGCTCGATAACAGATTTAATGAACTGCTTACGCGTCTTCCACCATCGGCTGCACCTGCCGCACCTCTGCAACAACTACTACTACTACTATACCTCCAGATCGCGAAACAGCCCTCCGCCGAGCGAGCCGTGTCCTTCTTCAGCCTGGCCAAACTGTTGGTGCTGCTGTTGATAGTTCTGCTGCTGATGCGGAGGGTGATTATGCGGGAGATTACGAGG

Zper:
>testassemblies4303(long)
CAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACACCCGTTTTCGGAATGGGTGACGTGCGGCAACAAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCTGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACACACGTTTTCGGAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAAAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCTGACACCCGTTTTCGGAATGGGTGACGTGCGG
>testassemblies2616
CAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACACCCGTTTTCGGAATGGGTGACGTGCGG

#Creating a bigger net: Grabbing all Cent assemblies from TRF finds

Unfortunately, the single monomer strategy is not sufficient to capture the diversity of
centromere repeats in all species.  Therefore, we have devised a new strategy of gathering
the diversity of repeats in the TRF assemblies.  We will make blast databases out of all of
the TRF finds, and blast the monomer from the most abundant polymer against them.  Those
with high blastn homology will be pulled out, and used as a larger mapping database.  All
of this is pursued on farm, /home/pbilinsk/huffwork/BLAST_TRF_Parsing/BiggerNet

makeblastdb -in Apludarepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Phylorepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Rimmarepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Sorghumrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Tanderrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Tdactrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Tperurepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Triturrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in Zperrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids

blastn -query centapluda.fa -evalue 1E-1 -outfmt 7 -db Apludarepeat_TRFfinds.fasta -task blastn -out DB_apluda
blastn -query centphylo.fa -evalue 1E-1 -outfmt 7 -db Phylorepeat_TRFfinds.fasta -task blastn -out DB_phylo
blastn -query centrimma.fa -evalue 1E-1 -outfmt 7 -db Rimmarepeat_TRFfinds.fasta -task blastn -out DB_rimmma
blastn -query centsorghum.fa -evalue 1E-1 -outfmt 7 -db Sorghumrepeat_TRFfinds.fasta -task blastn -out DB_sorgh
blastn -query centtander.fa -evalue 1E-1 -outfmt 7 -db Tanderrepeat_TRFfinds.fasta -task blastn -out DB_tander
blastn -query centtdact.fa -evalue 1E-1 -outfmt 7 -db Tdactrepeat_TRFfinds.fasta -task blastn -out DB_tdact
blastn -query centtperu.fa -evalue 1E-1 -outfmt 7 -db Tperurepeat_TRFfinds.fasta -task blastn -out DB_tperu
blastn -query centtritur.fa -evalue 1E-1 -outfmt 7 -db Triturrepeat_TRFfinds.fasta -task blastn -out DB_tritur
blastn -query centzper.fa -evalue 1E-1 -outfmt 7 -db Zperrepeat_TRFfinds.fasta -task blastn -out DB_zper

perl Blast_DBparser.pl DB_apluda | uniq | sed 's/test/>test/g' > bignetnames_apluda.txt
perl Blast_DBparser.pl DB_phylo | uniq | sed 's/test/>test/g' > bignetnames_phylo.txt
perl Blast_DBparser.pl DB_rimmma | uniq | sed 's/test/>test/g' > bignetnames_rimma.txt
perl Blast_DBparser.pl DB_sorgh | uniq | sed 's/test/>test/g' > bignetnames_sorghum.txt
perl Blast_DBparser.pl DB_tander | uniq | sed 's/test/>test/g' > bignetnames_tander.txt
perl Blast_DBparser.pl DB_tdact | uniq | sed 's/test/>test/g' > bignetnames_tdact.txt
perl Blast_DBparser.pl DB_tperu | uniq | sed 's/test/>test/g' > bignetnames_tperu.txt
perl Blast_DBparser.pl DB_tritur | uniq | sed 's/test/>test/g' > bignetnames_tritur.txt
perl Blast_DBparser.pl DB_zper | uniq | sed 's/test/>test/g' > bignetnames_zper.txt

perl Fetch_assembles.pl bignetnames_apluda.txt Apludarepeat_TRFfinds.fasta > Bignet_apluda.fa
perl Fetch_assembles.pl bignetnames_phylo.txt Phylorepeat_TRFfinds.fasta > Bignet_phylo.fa
perl Fetch_assembles.pl bignetnames_rimma.txt Rimmarepeat_TRFfinds.fasta > Bignet_rimma.fa
perl Fetch_assembles.pl bignetnames_sorghum.txt Sorghumrepeat_TRFfinds.fasta > Bignet_sorghum.fa
perl Fetch_assembles.pl bignetnames_tander.txt Tanderrepeat_TRFfinds.fasta > Bignet_tander.fa
perl Fetch_assembles.pl bignetnames_tdact.txt Tdactrepeat_TRFfinds.fasta > Bignet_tdact.fa
perl Fetch_assembles.pl bignetnames_tperu.txt Tperurepeat_TRFfinds.fasta > Bignet_tperu.fa
perl Fetch_assembles.pl bignetnames_tritur.txt Triturrepeat_TRFfinds.fasta > Bignet_tritur.fa
perl Fetch_assembles.pl bignetnames_zper.txt Zperrepeat_TRFfinds.fasta > Bignet_zper.fa

With the broader references made, I wanted to test what the % difference in mapping was when
using 100 sequences instead of the full 273 sequences.  This would give me an idea of how
badly biased the different-sized libraries could be in getting the diversity of centromere
repeat reads. I saw no major difference in the number of reads mapping to the full or half
size refs. Also, I wanted to test how well our de novo contigs were detecting centc compared
to our known centromere repeat library in maize.  I mapped against the single de novo monomer,
the bignet polymer, and the B73v2.centcseq reference, which is all 12,162 centc's in the v2
reference.  The monomer and the reference had approximately the same number of reads mapping
(close to 2000) while the polymer had 2.5x more (roughly 5000).  This suggests that our de
novo polymers are capturing a broader diversity of reads than are found in the reference.
Taken together, I believe the bignet polymers are a suitable, diverse reference to map against,
and, beyond a minimal read threshold, overall read depth does not have an appreciable impact
on the diversity captured in the bignet reference.

Make the cent and bignet references in the References folder, and map the library against
these two new references to get the broad and narrow abundance of the centromere repeat.

I think thats it.
	








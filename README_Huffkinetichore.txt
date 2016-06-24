Centromere Repeat Evolution in the Grasses

October 2013-February 2014
by: Paul Bilinski

to do things interactively: srun --pty -p bigmemh bash


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
panrepeat_out.unpadded.fasta

Round2
anepalRepeat_out.unpadded.fasta
hyphiRepeat_out.unpadded.fasta
isrugRepeat_out.unpadded.fasta
osatRepeat_out.unpadded.fasta
tflorRepeat_out.unpadded.fasta
tlaxRepeat_out.unpadded.fasta
udigRepeat_out.unpadded.fasta

Round3 - had to change manifest file for new version of mira: 
parameters = --hirep_good -NW:cnfs=no -NW:mrnl=200 -HS:mnr=no
tefRepeat_out.unpadded.fasta


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
./trf407b.macos64 panrepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h

Round 2
./trf407b.macos64 anepalRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 hyphiRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 isrugRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 osatRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 tflorRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 tlaxRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h
./trf407b.macos64 udigRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h

Round3
./trf407b.macos64 tefRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h


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
perl TRF_parser.pl panrepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > panrepeat_TRFfinds.fasta
#panicum has too few assemblies made.  not going to be a fair comparison.

Round2
perl TRF_parser.pl anepalRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > anepalrepeat_TRFfinds.fasta
perl TRF_parser.pl hyphiRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > hyphirepeat_TRFfinds.fasta
perl TRF_parser.pl isrugRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > isrugrepeat_TRFfinds.fasta
perl TRF_parser.pl osatRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > osatrepeat_TRFfinds.fasta
perl TRF_parser.pl tflorRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > tflorrepeat_TRFfinds.fasta
perl TRF_parser.pl tlaxRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > tlaxrepeat_TRFfinds.fasta
perl TRF_parser.pl udigRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > udigrepeat_TRFfinds.fasta

Round3
perl TRF_parser.pl tefRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > tefrepeat_TRFfinds.fasta

This will produce a file where each line is a tandem repeat and its name is >textassembly#

--- Abundance of non-Knob elements, what is in the top 10?

First, we want to understand what the most abundant of the tandem repeats is.  By far, the
most common element with be the knob repeat, and therefore we have to mask out this sequence.  
To do this, we use BLAST and the 34 annotated knob sequences from NCBI. Done on cluster,
in the directory: ~/huffwork/BLAST_TRF_Parsing/BiggerNet

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

Round2
blastn -query anepalrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_anepaltrf
blastn -query hyphirepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_hyphitrf
blastn -query isrugrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_isrugtrf
blastn -query osatrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_osattrf
blastn -query tflorrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_tflortrf
blastn -query tlaxrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_tlaxtrf
blastn -query udigrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_udigtrf

Round3
blastn -query tefrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_teftrf

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

Round2
perl Blast_DBparser.pl DB_knobs_vs_anepaltrf > assemblynamesinDB_knobs_vs_anepaltrf.txt
perl Blast_DBparser.pl DB_knobs_vs_hyphitrf > assemblynamesinDB_knobs_vs_hyphitrf.txt
perl Blast_DBparser.pl DB_knobs_vs_isrugtrf > assemblynamesinDB_knobs_vs_isrugtrf.txt
perl Blast_DBparser.pl DB_knobs_vs_osattrf > assemblynamesinDB_knobs_vs_osattrf.txt
perl Blast_DBparser.pl DB_knobs_vs_tflortrf > assemblynamesinDB_knobs_vs_tflortrf.txt
perl Blast_DBparser.pl DB_knobs_vs_tlaxtrf > assemblynamesinDB_knobs_vs_tlaxtrf.txt
perl Blast_DBparser.pl DB_knobs_vs_udigtrf > assemblynamesinDB_knobs_vs_udigtrf.txt

Round3
perl Blast_DBparser.pl DB_knobs_vs_teftrf > assemblynamesinDB_knobs_vs_teftrf.txt
#tef has no knob, cool

Notes: tdact has no hits (because its actually tritur).  sorghum has 1 hit.  The rest have a ton.  Those 2 were done by
hand for ReadyForMosaik files.  Second note, phylo has no knobs.

#ANEPAL and its unique knob.  Side project to filter out the unique knob from Anepal
#filtering done in biggernet, dont with ReadyForMosaik_anepal.txt so that it has no knobs

makeblastdb -in Anepal_knobs.txt -dbtype 'nucl' -parse_seqids
blastn -query ReadyForMosaik_anepal.txt -evalue 1E-1 -outfmt 7 -db Anepal_knobs.txt -task blastn -out DB_NEWknobs_vs_anepaltrf
Transfer files to my computer, the playdir directory, and run: 
perl Blast_DBparser.pl DB_NEWknobs_vs_anepaltrf > assemblynamesinDB_NEWknobs_vs_anepaltrf.txt
Open in textwrangler, process duplicates leaving 1.
load into R, find those that are not in the new knob.
perl Fetch_nonrepeatassemblies.pl Anepal_nonNEWknobassemblies.txt anepalrepeat_TRFfinds.fasta > test.txt
last one is messed up as in other files, change it over.
final product is: ReadyForMosaik_anepalnoNEWknob.txt

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

Round2
Notes: anepal and udig has a few, hand removed.  hyphi isrug osat have none.  tflor tlax
have a lot, removed via script.

grep ">" tflorrepeat_TRFfinds.fasta > tflor_assemblynames_forsubset.txt
grep ">" tlaxrepeat_TRFfinds.fasta > tlax_assemblynames_forsubset.txt

Round3
tef has no knobs, map away

Don't forget to get rid of the > in text wrangler.

The list is then imported into R, where I load 
the list of names of all of the test assemblies.  Use the short script:

	setwd("~/Documents/Projects/Huff_CentromereEvo/Git_Huff_CentEvo/playdir/")

	Orig <- read.csv("tflor_assemblynames_forsubset.txt", header=FALSE)
	remain <- read.csv("assemblynamesinDB_knobs_vs_tflortrf.txt", header=FALSE)
	notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
	write.csv(notknownrepeats,file = "nonknob_tflor.csv")

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

Round2
perl Fetch_nonrepeatassemblies.pl tflor_nonknobassemblynames.txt tflorrepeat_TRFfinds.fasta > test.txt
perl Fetch_nonrepeatassemblies.pl tlax_nonknobassemblynames.txt tlaxrepeat_TRFfinds.fasta > test.txt

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
#Round2: ran script on tflor and tlax, udig and anepal done by hand. rest no knobs.

Renamed the file since it is now ready for mosaik!  Moved onto the cluster, and mapped.

###Abundance of the De novo sequence 

---  Identifying the most common tandem repeat via mosaik

After renaming, moved to the cluster and made it into mosaik ready .dat.  Operations were
executed in the ~/huffwork/Mosaikmapping/ directory on the cluster.  Then build the 
sequences and submit the script via (all files must go to correct dir, and you have to
copy over the Submitfulltest.sh and the Mosaik executables, submit must also be altered
to contain bigmemh Q):

	./MosaikBuild -fr ReadyForMosaik_NAME.txt -oa ReadyForMosaik_NAME.dat
	sbatch Submitfulltest.sh ReadyForMosaik_nonknob.dat

The references go in the ~/huffwork/References directory, and make sure the fastq files are built for mosaik.
Once mapped, in each separate folder the alignments produce a *.loc file.  This will give a listing of each assembly and how many reads
mapped to it.  Do a head ___ > RanksNAME.csv  for each of the loc files.  Download them to
a local directory, and open them up in text wrangler.  Format this so you can open it as a
csv in excel.  Order them, find the top hit, and locate it in the TRF output file.  Find
the monomer version of the longer hit, those are shown below.

#anepal has some form of knob.  2 variants, a 185 and a 369 with no blast homology between them.
#anepal will have to be rerun after we filter out the knobs.
#Found high percentage of the genome mapping by lowering mhp to 1 so the run could finish

#isrug assemblies look unlike any of the others... tons of TRF finds in a single contig
#so again reduced the number of mhp to 1 to get a run to complete, saw 2% of the genome mapping
#took the top hit of the mhp1 run, it was a 185bp repeat with rough similarity to sorghum
#will use this to blast out all sequences, and see if i can still use the mhp of 100

Apluda: like sorghum, not
>testassemblies383(long)
GAAACTCATTTCGGCCTGTTTGGAGACTCAGATTTATCCCAGTGCAACATAGGTGCACGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGTGGAAACTAGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCAGTGCAACATAGGTGCCCGGTTTGCATCGAATGTACCATCAGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGTGGAAACTCGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCAGTGCAACATAGGTGCCCGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGTGGAAACTCGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCGGTGCAACATAGGTGCCCGATTTGCGTCCAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTACTAACGCGGGTCCAGTGGAAACTCGTTTTGGCCTGTTTGGAGACTCAGATTTATCCCGGTGCAACATAGGTGCCCGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAGT
>testassemblies384(short, 138bp)
TGGAAACTCGTTTCGGCCTGTTTGGAGACTCAGATTTATCCCGGTGCAACATAGGTGCCCGGTTTGCGTCGAATGTACCATAGGCTTGGGACTCACGCTGGGCACACCCTATGGTACTCCGTAGTAACGCGGGTCCAG

**Phylo: blasts to bamboo repetitive region, not in melters, not annotated as anything
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

Zper: maize centc
>testassemblies4303(long)
CAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACACCCGTTTTCGGAATGGGTGACGTGCGGCAACAAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCTGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACACACGTTTTCGGAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAAAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCTGACACCCGTTTTCGGAATGGGTGACGTGCGG
>testassemblies2616
CAACGAAATTGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCGTGCTTTGTATGCACCCCGACACCCGTTTTCGGAATGGGTGACGTGCGG

ROUND2:

anepal: (has a knob like thing, 185bp 369bp)
>testassemblies1321(long, 504)
CATCAAAACGAGATTGCACACGACCCACATCATCTGGAAGTACAATCGGATGCATGCAAACCGAGCATTAAGTTGGCACCGAGATTAACACAGTCACCAAACGCATAAAAACGAGATTGCACACGACCCACATCATCTGGGAGTACAACCGAATGCTTACAAACAGACCTCTAACTTAGCACCGAGATTAACACAGTCTCCAAACGCATCAAAACAAGATTGCACACGACCCACGTCATCTAGAAGTACAATCGGATGCATGCAAACCGAGCTCTAAGTTTCCACCGAGATTAACACACTCTCCAAACGCATCAAAATGAGATTGCACACGACCCACGTCATCTGGAAGTACAATCGGATGCATGCAAACAGAGCTCTAAGTTTGCACCGAGATTAACACAGTCTCCAAACGCATCAAAAGAGATTGCACACGACTCACGTCATTTGGGAGTACATCGGTGCATGCAACAGAGCTCTAAGTTTCACCGAGATTAACATAGTCG
>testassemblies924(short, 103bp repeat)
AAAACGAGATTGCACACGACCCACGTCATCTGGGAGTACAATCGGATGCATGCAAACAGAGCTCTAAGTTTGCACCGAGATTAACACAGTCTCCAAACGCATC

#calling of the anepal repeat comes from homology to saccharum/sorghum.  chunk is repeated for 79bp, and there are 24bp till the next repeated section)

isrug:
>testassemblies2789(137bp)
ACCTTGGAGTATTGTCGGGTGTGCTAAAATTGCTCCCGCTCCAACGGGAGCTTCAGAACAACCCGTGCACATATTTTGTGCCGAAATTCACACGAGTCTCTAAATGAACCGAAACGAGCTTCCACCTCACCCACGTC
short, twice as abundant as next up which is TE

tlax: maize like
>testassemblies176(long)
TGAGTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATGGCGTGAAACCACCCCAAACATGAGTTTTGTACCTAATTAGTGGATTGGGCATGTTCGTTGCAAAAAAAAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCTGCAACGAAATGGCGCAAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTCCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTCCGCAATGGGTGACGTGCAACAACGAAATGGCGCGAAACCACCCCAAACATGAGTTTTGGAACTAAAGTAGTGGATTGGGCATGTTCGTTACGAAAAACGAAGAAATGGTTCCGGTGGCAAAAACTCATGCTTTGTATGCATCCCGATACCCGTTTTCGGAATAGGTGACGTGCGGCAACGAAATGGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTCGAATGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACCCATGCCTTGCATGCACACCGATACCCGTTTCCGAAATGGGTGACGCGCGGCAACGAAATGGCACGAAACCACCCCAAACA
>testassemblies175(short, 156bp, maize like)
TGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGGAATGGGTGACGTGCGGCAACGAAATGGCGCGAAACCACCCCAAACA

**udig:
>testassemblies613(long)
AAACGTACCCTGGCCTCTTCACAAGATTGCTAGTCATCAATTTGAGGCACTGGCCCTGAGAGAAATTGGCAACACAGTGCCACCATGGAGGGTTTTGTCCACCCAAGCCTGTTTAGGCCATTTTGGGACTCAGTTTTTGCTTCGTGCGATACTATGTTTCTCAAAATGGTTTTTAGTTATAGCAAAACATACCCTGGCCTCTTCACAAATAGCTAGTCATCAATTTGAGGCACTGGCCCTGAGAGAAATTGGCAACACGGTGCCGCATGGAGGGTTTTGTCCACCCAAGCCTGTTTAGGCTGTTTTGGGACTCAGTTTTTGCTTCGTGCGATACTATGTTTCTCAAAATAGTTTTTGGTTATAGAGAAACATACCCTGGCCTCTTCACATGATAGCTAGTCATCAATTTGAGGCACTGGCCCTGAGAGAAATTGGCAACACGGTGCCGCCATGGAGGGTTTTGTCCACCCAAGCCTGTTTAGGCCGTTTTGGGACTCAGTTTTTGTTTAGGCCATTTTGGGACTCAGTTTTTGCTTCGTGCAATACTATGTGTCTCAAAATAGTTTTTGGTTATAGCAAATAATACCCTAGCCTCTTCACACGATAGCTAGTCATCAATTTGAGGCACTAGCCCTAAGAGAAATTGCCAACACGAAGCCGCCATGGAGGGTTTTGTCAACCCAAGCCTGTTTAGGCCGTTTTGGGACTCAGTTTTTGCTTCGTGCGATACTATGTTTCCCAAAATGGTTTTTGGTTATAGCAAATAATACCCTGGCCTCTTCACACAATAGCTAGTCATCAATTTGAGGCACTAGCCTAGAGAGAAATTGGCAACACGGTGCCGCCATGGAGGGTTTTGTCCACCCAAGCCTATTTAGGCCATTTTGGGACGCAGTTTTTGATTCGTGCGATACTATGTTTCTCAAAAAAGTTTTTGGTTATAGCGAAAAGTACCCTAGCCACATCACACGATAGAAACTCATCAATTTGAGGCACTGGCCCTGAGAGAAATTGGCAACACGGTGCCGCCATGGAGGGTTTTGTCCACCCAAGCCTGTTTAGGCCATTTTGGGACTCAGTTTTTGCTTTGTGCGATACTATGTTTCTCAAAATGGTTTTTGGGTATAGCG
>testassemblies605(short, 183bp, no homology)
GTTATAGCAAAACATACCCTGGCCTCTTCACATGATAGCTAGTCATCAATTTGAGGCACTGGCCCTGAGAGAAATTGGCAACACGGTGCCGCCACGAAGGGTTTTGTCCACCCAAGCCTGTTTTGGCTATTTTGGGACTCAGTTTTTGCTTCGTGCGATACTATGTTTCTCAAATGGTTTTTA

**hyphi: brand new, 157bp long
>testassemblies521
GATTTTTCGTGCGAAAGTGGCGTGCACGGGCATCCCCGTCGAGTTTTGCACGTTTTTGCTCCGAACGAGATCCGAAAGCGCGAAACGCGCCAAACACGTGTCCCCCAACCAAACAGGGTTGGATGTGTGCGTTCGTTGCGAAAAATTGCGCCGCACGATTTTTCGTGCGGAAGCGGCGTGCACGGGCATCCCCGTCAAGTTTCGCACGTTCTTGCTCAGAACGAGATCCGAAAGTCGCGAAACGCGCCAAACACGTCTCCCCCAACCAAACGGGGTTGGATGTGTGCGTTCATTGCGAAAAATTGCACCGCACGATTTTCGTGCGACAATGGCGTGCACGGGCATCCCCGTCGAGTTTCGCACGTTCTTGCTCCGATCGAGATCCGAAAGACGCGAAACGCGCCAAACACGTGTCCCCAACCAAACAGGGTTGGATGTGTGCGTTCGTTGCGAAAAATTGCGCCGCACGATTTTTCGTGCGGAAGTGGCATGCACAGGCATCCCCGTCGAGTTTCGCACGTTCTTGCTCCGAACGAGATCCGAAAGACGAGAAACGCGCCAAACACGTCTCCCCCAACCAAACGGAGTTGGATGTGTGCGTTCGTTGCGAAAAATTGCACCGCACGACTTTTCGTGCGAATGTGGCGTGCACGGGCATCCCCGTCGAGTTTCGCACGTTTTTGCTCCGAACGAGATCCGAAAGATGCGAAACGCGCCAAACACATGTCCCCCAACCAAACGGGGTTGGACGTGTGCATTCGTTGCGAAAAATTGCACCGCACAATTTTTCGTGCGACAGTGGCGTACACGGGCATCCCCGTCGAGTTTCGCACGTTCTTGCTCCGAACGACATCCGAAAGACGCGAAACACGCCAAACACGTGTCCCCCAACCAAACGGAGCTGGATGTGTGCGTTCGTTGCGAAAAATTACACCGCAC
>testassemblies520(short 157bp no homology to anything)
GATTTTTCGTGCGGAAGTGGCGTGCACGGGCATCCCCGTCGAGTTTCGCACGTTCTTGCTCCGAACGAGATCCGAAAGTCGCGAAACGCGCCAAACACGTGTCCCCCAACCAAACGGGGTTGGATGTGTGCGTTCGTTGCGAAAAATTGCGCCGCAC

tflor: just like maize
>testassemblies2478(long)
ATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATGGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAATAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATTGCGCGAAACCACCCCAAACATGAGTATTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATGGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATGGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGGCAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATGGCGCAAAAACACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCAGTGGCAAAAACTCATGCCTTGTTATGCACCACGATACCCATTTTCGCAATGGGTGACGTGCGGCAACGAAATGACGCAAAACCACCCCAAACATGAGTTTTGGACCTAAAGCAGTGGATTGGGCATGTTCGTTGCGAAAAAAAGAAATGGTTCCGGTGGCAAAAACTCATGCATACA
>testassemblies2477(short 156bp)
CAAAAACTCATGCCTTGTATGCACCCCGATACCCGTTTTCGCAATGGGTGACGTGCGGCAACGAAATGGCGCGAAACCACCCCAAACATGAGTTTTGGACCTAAAGTAGTGGATTGGGCATGTTCGTTGCGAAAAAAGAAGAAATGGTTCCGGTGG

osat: all 10 are RCS2 like, but not quite as long.  have homology over entire length tho.
>testassemblies8(long)
CTCACTTCGTGATTCGCGCGGCGAACTTTTGTCAATTAATGCCAATATTGGCACACGAGGGTGCGATGTTTTTGACCGGAATCAAAAAGTTCAAAAAAACCAAAACATGATTTTTGGACATATTGGAGTGTATTGGGTGCGTTCGTGGAAAAAACTCACTTCGTGATTCGCGCGGCGAACTTTTGTCAATTGATGCCAATATTGGCACACAGGGTGCGATGTTTTTGACCGGAATCAAAAAGTTCGAAAAAAAACCAAAACATGATTTTTGGACATATTGGAGTGTATTGGGTGCGTTCGTGGCAAAAACTCACTTCGTGATTCGCGCGGCGAACTTTTGTCAATTGATGCCAATATTGGCACACGAGGGTGCGATGTTTTTGTCTGGAATCAAAAAGTTCAAGAAAAACGAAACATGATTTTTGGATATATTGGACTGTATTGGGTGCGTTCGTGGCAAAAACTCACTTCGTGATTCGCGCGGTGAACTTTTGTCAATTGATGCCAATATTGGCACACGAGGGTGCGATGTTTTTGACCGGAATCAAAAAGTTCAAAAAAAACCAAAACATGATTTTTGGACATATTGGAGTGTATTGGGTGCGTTCGTGGCAAAAACTCACTTCGTGATTCGCGCGGCGAACTTTTGTCAATTATGCCAATATTGGCACACGAGGTGTGATGTTTTTGACCGGAATCAAAAAGTTAAAAAAAAACCAAAACATGATTTTTGGACATATTGGAGTGTATTGGGTGCGTTCGTGGCAAAAA
>testassemblies6(short, 154bp)
GGTGCGATGTTTTTGACCGGAATCAAAAAGTTCAAAAAAAACAAAACATGATTTTTGGACATATTGGAGTGTATTGGGTGCGTTCGTGGCAAAAACTCACTTCGTGATTCGCGCGGCGAACTTTTGTCAATTGATGCCAATATTGGCACACGAG

Round3
tef: brand new, genome paper does not have any centromere repeat annotated as far as i could tell
>testassemblies525(long, no real homology to anything, multiple shorter repeats)
GTAGAGCTCCGAATCAGCCTGGTTTGACTCGTTGGCAACTAAACGCAAGTTTTTGAGTGTTTTCCCATGGAAACCATCCAATTACATCAAAAACATCATAACAAGTGAAAAAGAACGATACTCGGGTCGTTTTGACCGAAACAACACGGCATTCACCCAAACGGCTCCCGGAGAGATCCGAATCAGCCCGGTTTTCACAAGTTTGCAACTAAACACAAGTTTTTGAGTGTTTTCCCATGGGAACCATCCAATTACATCAAAAACATCATAACAAGTGAAAAAGAACGATTCACTCGGGACGTTTTGACCGAACAACGATACCAACGGCTCCCGTAGAGCTCCGAATCAGCCCGGTTTTGACTCGTTTGCAACTAAACACAAGTTTTTGAGTGTTTTCCCATGGGAACCATCCAATTACATCAAAAACATCATAACAAGTGAAAAAAATGATACTCGGGTCGTTTTGACCGAAACAACACGATATTCACCCAAACGGCTCCCGTAGAGCTCCAAATCAGCCCGGTTTGACTCATCAAAAAAATCATTAGTAAAAAGATGATACTCGGTCGTCTTGACCGAAGCAACACGGAATTCACAAAAACAGCTCAAGTAGAGCTCCGAATCAGCCCGGTTTGACTCGTTTGCAACTAAACACAAGTTTTGAGTGTTTTCCCATGGGAAACAATCAAATTAAATCAAAAACATCATAACAAGTGAAAAAGAATGATACTCGGGTCATTTTGACCAAAACAACACGGTATTCACCCAAACGGCTCCCTAGAGCTCCGAATCAGCCCGGTTTGACTCCTTTGCAACTAAACACAAGTTTTTGAGTGTGCTCACATAGGAACCATCCAACTAAGTCAAAAACATCATAAAAAGTGAAAAAGAATGATTCACTCGGGTCGTTTTGACCGAAACAACACGGTATTCACCAAAACGGCTCCCGTAGAGCTCCGAACCAGCCCGGTTTGACTCATTTGCAACTAAACACAAGTTTTTGAGTGTTTTCCCATGGGAACCATCCAATTACATCAAAAAAAGCATAACAAGTGAAAAAGAATGATACTCGGGTCGTTTTGACCAAAACAACACGGTATTCACCCAAACGGCTCCCGTAGAGCTCCGAATCAGCCCGGTTTGACTCCTTTGCAAGTAAACACTTGCTTTTGAGTGTTTTCCCATGAGAACCATCCAATTACATCAAAAAAATCATAACAAGTGAAAAAGAATGATACTCGGGTCGTTTTGACCAAAACAACACGGTATTCACCTGAGTGTTTTCCTATGGGAACCATCCAATTACATCAAAAACATCATAACAAGCTA
>testassemblies517(short, 168bp consensus, was most abundant repeat)
GGTCGTTTTGACCGAAACAACACGGTATTCACCCAAACGGCTCCCGTAGAGCTCCGAATCAGCCCGGTTTGACTCGTTTGCAAGTAAACAAGTTTTTGAGTGTTTTCCCATGGGAACCATCCAATTACATCAAAAACATCATAACAAGTGAAAAAGAATGATACTCG

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

Round2:
makeblastdb -in anepalrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in tlaxrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in tflorrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in osatrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in isrugrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in hyphirepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
makeblastdb -in udigrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids

blastn -query centanepal.fa -evalue 1E-1 -outfmt 7 -db anepalrepeat_TRFfinds.fasta -task blastn -out DB_anepal
blastn -query centisrug.fa -evalue 1E-1 -outfmt 7 -db isrugrepeat_TRFfinds.fasta -task blastn -out DB_isrug
blastn -query centtlax.fa -evalue 1E-1 -outfmt 7 -db tlaxrepeat_TRFfinds.fasta -task blastn -out DB_tlax
blastn -query centudig.fa -evalue 1E-1 -outfmt 7 -db udigrepeat_TRFfinds.fasta -task blastn -out DB_udig
blastn -query centhyphi.fa -evalue 1E-1 -outfmt 7 -db hyphirepeat_TRFfinds.fasta -task blastn -out DB_hyphi
blastn -query centtflor.fa -evalue 1E-1 -outfmt 7 -db tflorrepeat_TRFfinds.fasta -task blastn -out DB_tflor
blastn -query centosat.fa -evalue 1E-1 -outfmt 7 -db osatrepeat_TRFfinds.fasta -task blastn -out DB_osat

perl Blast_DBparser.pl DB_anepal | uniq | sed 's/test/>test/g' > bignetnames_anepal.txt
perl Blast_DBparser.pl DB_isrug | uniq | sed 's/test/>test/g' > bignetnames_isrug.txt
perl Blast_DBparser.pl DB_tlax | uniq | sed 's/test/>test/g' > bignetnames_tlax.txt
perl Blast_DBparser.pl DB_udig | uniq | sed 's/test/>test/g' > bignetnames_udig.txt
perl Blast_DBparser.pl DB_hyphi | uniq | sed 's/test/>test/g' > bignetnames_hyphi.txt
perl Blast_DBparser.pl DB_tflor | uniq | sed 's/test/>test/g' > bignetnames_tflor.txt
perl Blast_DBparser.pl DB_osat | uniq | sed 's/test/>test/g' > bignetnames_osat.txt

perl Fetch_assembles.pl bignetnames_anepal.txt anepalrepeat_TRFfinds.fasta > Bignet_anepal.fa
perl Fetch_assembles.pl bignetnames_isrug.txt isrugrepeat_TRFfinds.fasta > Bignet_isrug.fa
perl Fetch_assembles.pl bignetnames_tlax.txt tlaxrepeat_TRFfinds.fasta > Bignet_tlax.fa
perl Fetch_assembles.pl bignetnames_udig.txt udigrepeat_TRFfinds.fasta > Bignet_udig.fa
perl Fetch_assembles.pl bignetnames_hyphi.txt hyphirepeat_TRFfinds.fasta > Bignet_hyphi.fa
perl Fetch_assembles.pl bignetnames_tflor.txt tflorrepeat_TRFfinds.fasta > Bignet_tflor.fa
perl Fetch_assembles.pl bignetnames_osat.txt osatrepeat_TRFfinds.fasta > Bignet_osat.fa

Round3
makeblastdb -in tefrepeat_TRFfinds.fasta -dbtype 'nucl' -parse_seqids
blastn -query centtef.fa -evalue 1E-1 -outfmt 7 -db tefrepeat_TRFfinds.fasta -task blastn -out DB_tef
perl Blast_DBparser.pl DB_tef | uniq | sed 's/test/>test/g' > bignetnames_tef.txt
perl Fetch_assembles.pl bignetnames_tef.txt tefrepeat_TRFfinds.fasta > Bignet_tef.fa
#Checking position in the tef genome:  blast, has 416 hits across all scaffolds, mostly to
position on scaffolds.  Not super meaningful.  However, we take the bigger scaffolds, those
greater than 500,000bp, not a lot of meaningful hits.  The larger scaffolds can be found in
the file Etbiggerscaffold.fa

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
	








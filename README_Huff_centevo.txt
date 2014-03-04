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
in using either rice or maize or sorghum (see CRR CRS CRM testing in PB directory).  To
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

--- Using Tandem repeat finder on the assembly file from mira and filter

In the _results folder, you will find the outputs of the assembly.  One of them will be
the unpadded file.  example:

RimmarepeatASS_out.unpadded.fasta

Download this file from the cluster, and run trf on it locally.

./trf407b.macos64 RimmarepeatASS_out.unpadded.fasta 2 7 7 80 10 50 2000 -h

Next, take the .dat produced by TRF and filter it based on the length requirements we see
in Melters et al 2013.  They found that the shortest tandem repeat in a plant was 40bp
long, in ricin.  The script TRF_parser.pl will extract all sequences whose tandem repeat
length is more than 40bp.  

perl TRF_parser.pl RimmarepeatASS_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Rimmarepeat_TRFfinds.fasta

This will produce a file where each line is a tandem repeat and its name is >textassembly#

--- Abundance of non-Knob elements, what is in the top 10?

First, we want to understand what the most abundant of the tandem repeats is.  By far, the
most common element with be the knob repeat, and therefore we have to mask out this sequence.  
To do this, we use BLAST and the 34 annotated knob sequences from NCBI.  

makeblastdb -in 34KnobsonNCBI_renamed.txt -dbtype 'nucl' -parse_seqids

With the DB made, I then blast our tandem repeats against the database:

blastn -query Rimmarepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_rimmatrf

Transfer files to my computer, the playdir directory, and run: 

perl Blast_DBparser.pl DB_knobs_vs_rimmatrf > assemblynamesinDB_knobs_vs_rimmatrf.txt

I print out the names of each sequence that has a hit longer than 30bp to one of the 
repeats.  This is printed out, duplicates are removed via textwrangler (process duplicate
lines, delete duplicate lines, print to new file named assemblynamesinDB_GRAMvsRIMMA_nodupl.txt).  ###DID NOT RENAME
The original names of all of the assemblies are taken using 

	grep ">" Rimmarepeat_TRFfinds.fasta > Original_assemblynames_forsubset.txt

Don't forget to get rid of the > in text wrangler.

The list is then imported into R, where I load 
the list of names of all of the test assemblies.  Use the short script:

	Orig <- read.csv("Original_assemblynames_forsubset.txt", header=FALSE)
	remain <- read.csv("assemblynamesinDB_GRAMvsRIMMA2_nodupl.txt", header=FALSE)
	notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
	write.csv(notknownrepeats,file = "nonrepeatsknob.csv")

This will create a csv file with additional columns and such.  I just created a new file
from these names by copy pasting all of them into nonrepeatIDtags.txt, and adding the > 
with replace of testassemblies for >testassemblies.  Now, we need to use those names that
dont have a hit against all the repeats to get the assemblies are unknown.  We have the
script Fetch_nonrepeatassemblies.pl.

	perl Fetch_nonrepeatassemblies.pl nonknobIDtags.txt Rimmarepeat_TRFfinds.fasta > test.txt

Then, test whether those names are the same in each of the files.

	grep ">" test.txt > compareafter
	grep ">" nonknobIDtags.txt > comparebefore
	diff comparebefore compareafter 

In test.txt, I found that the very last entry was wrong, and replaced it... Will have to do this for
each of the individuals.  I replaced:

	>testassemblies0
	CTGGCTCGGGCCGATTCCAGCGTAAACCGTGAGCTAAAACAGCGTAAATGAGTATAGAAATTTAGCGTAAATCTTATATCTGTTTTGTAACAGCAAATGAGGCCTAAAATTACGGCGTGAAAATTGTGTGCAGCCGATCGTGCACGGGTC

	>testassemblies5652
	AAATATGAAAATACAATTTAAATGAATTAATATTTTAATTTTTAACAG

	mv test.txt ReadyForMosaik_nonknobs.txt

Renamed the file since it is now ready for mosaik!  Moved onto the cluster, and mapped.

###Abundance of the De novo sequence 

---  Identifying the most common tandem repeat via mosaik

After renaming, moved to the cluster and made it into mosaik ready .dat.  Operations were
executed in the ~/huffwork/Mosaikmapping/NonreapeatRimma/ directory on the cluster.  Then
build the sequences and submit the script via (all files must go to correct dir):

	./MosaikBuild -fr ReadyForMosaik_nonknob.txt -oa ReadyForMosaik_nonknob.dat
	sbatch Submitfulltest.sh ReadyForMosaik_nonknob.dat


	sbatch Submitfulltest.sh ReadyForMosaik_nonknobs.dat
	in NonknobRimma

The top 10 hits after knobs are removed are:

17389	testassemblies2588 centc, 3 to 5 monomers
16708	testassemblies3069 centc, 4x
11442	testassemblies4509 centc, 3x
10081	testassemblies1231 centc, 2x
10037	testassemblies4883 centc, 2x
9765	testassemblies3915 centc, 2x
9673	testassemblies3288 centc, 1 to 2x
9284	testassemblies1267 centc 3x
8665	testassemblies2997 centc 3x
8315	testassemblies4508 centc 2x

But they have different lengths, and are imperfect assemblies of CentC.  Always, the longest
assembly will have the most hits because its the longest.

!!! Go back to the original sequence testassemblies2588 and see how long the assembly is
and how that compares to the TRF find from it.

### Re-do Centromere Repeat analysis for Sorghum

Submit script via sbatch Submit_assembly.sh


	project = SorghumRepeat 
	job = genome,denovo,accurate
	parameters = -highlyrepetitive -NW:cnfs=no -NW:mrnl=200 -HS:mnr=no

	readgroup = Rimma
	data = /home/pbilinsk/huffwork/Data/FWD/MH_Sorghum_bicolor_sp_bicolor.fastq
	technology = solexa


./trf407b.macos64 SorghumRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h

perl TRF_parser.pl SorghumRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Sorghumrepeat_TRFfinds.fasta

>sorghum_centromere_repeat_melters
GAGATAGTGTTAATCTTGATGCAAGATAGGTGCACGGTTTGCACGGAACGCACCATAGGCTAAGAAACCATTTTGGACGCACCCGATGGAACTCCTAGATGAAGTGTGTCAAATGGAAGCTCGGTTCGGTCTGTTTG 

blastn -query Sorghumrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_sorghumtrf

perl Blast_DBparser.pl DB_knobs_vs_sorghumtrf > assemblynamesinDB_knobs_vs_sorghumtrf.txt

Nothing that hits knob in sorghum.

mine -    Sorghum bicolor centromere-associated repetitive DNA element CEN38	 717	 2950	100%	0.0	90%	AF137608.1
melters - Sorghum bicolor centromere-associated repetitive DNA element CEN38	 167	 894	100%	3e-38	88%	AF137608.1

testassemblies13	 163	 614	100%	4e-45	87%	 57943


4367	testassemblies13 cent associated
2740	testassemblies151 same cent assoc
1887	testassemblies29 3rd is cent
1778	testassemblies12 3rd is cent
1737	testassemblies101 3rd is cent
1688	testassemblies53
1658	testassemblies132
1638	testassemblies51
1621	testassemblies98
1619	testassemblies16



Again, tripsacum dactyl

  project = TdactRepeat 
        job = genome,denovo,accurate
        parameters = -highlyrepetitive -NW:cnfs=no -NW:mrnl=200 -HS:mnr=no

        readgroup = TripD_mergefwdrev
        data = /home/pbilinsk/huffwork/Data/FWD/MH_Tripsacum_dactyloides.fastq
        technology = solexa
        
sbatch Submit_assembly.sh

./trf407b.macos64 TdactRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h

perl TRF_parser.pl TdactRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Tdactrepeat_TRFfinds.fasta

blastn -query Tdactrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_tdacttrf

Short homology to knobs, nothing over 20bp.  so nothing thrown out

./MosaikBuild -fr ReadyForMosaik_tdact.txt -oa ReadyForMosaik_tdact.dat
sbatch Submitfulltest.sh ReadyForMosaik_tdact.dat


Again, for tripsacum andersonii

./trf407b.macos64 TripAnderRepeat_out.unpadded.fasta 2 7 7 80 10 50 2000 -h

perl TRF_parser.pl TripAnderRepeat_out.unpadded.fasta.2.7.7.80.10.50.2000.dat > Tanderrepeat_TRFfinds.fasta

blastn -query Tanderrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_knobs_vs_tandertrf

perl Blast_DBparser.pl DB_knobs_vs_tandertrf > assemblynamesinDB_knobs_vs_tandertrf.txt

	grep ">" Tanderrepeat_TRFfinds.fasta > Tander_assemblynames_forsubset.txt

Don't forget to get rid of the > in text wrangler.

The list is then imported into R, where I load 
the list of names of all of the test assemblies.  Use the short script:

	Orig <- read.csv("Tander_assemblynames_forsubset.txt", header=FALSE)
	remain <- read.csv("assemblynamesinDB_knobs_vs_tandertrf.txt", header=FALSE)
	notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
	write.csv(notknownrepeats,file = "Tandernonrepeat.csv")

This will create a csv file with additional columns and such.  I just created a new file
from these names by copy pasting all of them into TandernonrepeatID.csv, and adding the > 
with replace of testassemblies for >testassemblies.  Now, we need to use those names that
dont have a hit against all the repeats to get the assemblies are unknown.  We have the
script Fetch_nonrepeatassemblies.pl.

	perl Fetch_nonrepeatassemblies.pl TandernonrepeatID.csv Tanderrepeat_TRFfinds.fasta > test.txt

Then, test whether those names are the same in each of the files.

	grep ">" test.txt > compareafter
	grep ">" TandernonrepeatID.csv > comparebefore
	diff comparebefore compareafter 

./MosaikBuild -fr ReadyForMosaik_Tander.txt -oa ReadyForMosaik_Tander.dat
	sbatch Submitfulltest.sh ReadyForMosaik_nonknob.dat



##BREAK


--- Removal of any tandem repeats that hit known Gramineae Repeats via BLAST

We want to remove any of the tandem repeats that BLAST to known repetitive elements.  To
get a reference for this, I downloaded all of the Gramineae repeats that are TEs, Telomere
, rDNA, 45SrDNA, 5SrDNA (from http://plantrepeats.plantbiology.msu.edu/downloads.html).
To this, I added all 34 knob sequences on NCBI.  The final file is located in
the folder Gramineae_repeats.  I built those sequences into a blast database with (in the
BLAST_TRF_Parsing folder):

makeblastdb -in 34KnobsonNCBI_renamed.txt -dbtype 'nucl' -parse_seqids

With the DB made, I then blast our tandem repeats against the database:

blastn -query TIGR_Gramineae_Repeats_nocent.fasta -evalue 1E-1 -outfmt 7 -db Rimmarepeat_TRFfinds.fasta -task blastn -out DB_GRAMvsRIMMAUGH

So when doing this, tons of knobs were left with megablast.  Adding -task blastn, so we can
get the knobs.  
This is a blastn run on the cluster that uses megablast, meaning it is somewhat strict. See
above for amendment to this.

Transfer files to my computer, the playdir directory, and run: 

perl Blast_DBparser.pl DB_GRAMvsRIMMA3 > assemblynamesinDB_GRAMvsRIMMA3.txt

I print out the names of each sequence that has a hit longer than 30bp to one of the 
repeats.  This is printed out, duplicates are removed via textwrangler (process duplicate
lines, delete duplicate lines, print to new file named assemblynamesinDB_GRAMvsRIMMA_nodupl.txt).
The original names of all of the assemblies are taken using 

	grep ">" Rimmarepeat_TRFfinds.fasta > Original_assemblynames_forsubset.txt

Don't forget to get rid of the > in text wrangler.

The list is then imported into R, where I load 
the list of names of all of the test assemblies.  Use the short script:

	Orig <- read.csv("Original_assemblynames_forsubset.txt", header=FALSE)
	remain <- read.csv("assemblynamesinDB_GRAMvsRIMMA2_nodupl.txt", header=FALSE)
	notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
	write.csv(notknownrepeats,file = "nonrepeats.csv")

This will create a csv file with additional columns and such.  I just created a new file
from these names by copy pasting all of them into nonrepeatIDtags.txt, and adding the > 
with replace of testassemblies for >testassemblies.  Now, we need to use those names that
dont have a hit against all the repeats to get the assemblies are unknown.  We have the
script Fetch_nonrepeatassemblies.pl.

	perl Fetch_nonrepeatassemblies.pl nonrepeatIDtags3.txt Rimmarepeat_TRFfinds.fasta > test.txt

Then, test whether those names are the same in each of the files.

	grep ">" test.txt > compareafter
	grep ">" nonrepeatIDtags3.txt > comparebefore
	diff comparebefore compareafter 

In test.txt, I found that the very last entry was wrong, and replaced it... Will have to do this for
each of the individuals.  I replaced:

	>testassemblies0
	CTGGCTCGGGCCGATTCCAGCGTAAACCGTGAGCTAAAACAGCGTAAATGAGTATAGAAATTTAGCGTAAATCTTATATCTGTTTTGTAACAGCAAATGAGGCCTAAAATTACGGCGTGAAAATTGTGTGCAGCCGATCGTGCACGGGTC

	>testassemblies5652
	AAATATGAAAATACAATTTAAATGAATTAATATTTTAATTTTTAACAG
	
This was the same for the nonrepeatIDtags2.txt iteration, the blastn go at it.

	mv test.txt ReadyForMosaik_nonrepeatassemblies3.txt
	
Renamed the file since it is now ready for mosaik!  Moved onto the cluster.

--- Does the library have CentC in it?  Using BLAST.

Before we start doing some mapping, just out of curiousity, we want to identify the 
assemblies that have the highest similarity to CentC, if there are any.  I grabbed all
12,162 CentC's in the maize refgen2, and made a blast reference from it.

makeblastdb -in test_allcentc.fasta -dbtype 'nucl' -parse_seqids

Then, blast all of our sequences that made it through the repeat filter (alternatively,
could use all of the sequences we have from the TRF finder)

blastn -query ReadyForMosaik_nonrepeatassemblies3.txt -evalue 1E-1 -outfmt 7 -db test_allcentc.fasta -num_threads 4 -out DB_NONREPEAT3_centc
blastn -query Tdactrepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db test_allcentc.fasta -num_threads 4 -out DB_TDACT_centc

Note that we have a ton of CentC before and after filtering.  CentC survives!

###Abundance of the De novo sequence 

---  Identifying the most common tandem repeat via mosaik

After renaming, moved to the cluster and made it into mosaik ready .dat.  Operations were
executed in the ~/huffwork/Mosaikmapping/NonreapeatRimma/ directory on the cluster.  Then
build the sequences and submit the script via (all files must go to correct dir):

	./MosaikBuild -fr ReadyForMosaik_nonrepeatassemblies.txt -oa ReadyForMosaik_nonrepeatassemblies.dat
	sbatch Submitfulltest.sh ReadyForMosaik_nonrepeatassemblies.dat

This process takes 3.5 days, so allow running time.  Move to folder:

cd /home/pbilinsk/huffwork/Mosaikmapping/Tests/NonreapeatRimma

The entire output of the will be in the folder.  The file you are interested in is .dat.loc
This will have the STDOUT for the last step of mosaik, which prints out how many reads map
to each of the reference sequences.  Figure out how many lines of repeats there are, print
to a csv file using head.  This way we can sort the read ID's by most common.  File name:

/home/pbilinsk/huffwork/Mosaikmapping/NonreapeatRimma/NonrepeatOutput.loc.txt

---
 
 
blastn -query knob.fasta -evalue 1E-1 -outfmt 7 -db trfseq.fasta -task blastn -out dbtrf

blastn -query trfseq.fasta -evalue 1E-1 -outfmt 7 -db knob.fasta -task blastn -out dbknob

blastn -query TIGR_Gramineae_Repeats_nocent.fasta -evalue 1E-1 -outfmt 7 -db Rimmarepeat_TRFfinds.fasta -task blastn -out DB_GRAMvsRIMMAUGH

blastn -query 34KnobsonNCBI_renamed.txt -evalue 100 -outfmt 7 -db Rimmarepeat_TRFfinds.fasta -task blastn -out DB_Knobsvstrfrimma


blastn -query 34KnobsonNCBI_renamed.txt -evalue 1E-1 -outfmt 7 -db trfseq.fasta -task blastn -out db34s



blastn -query Rimmarepeat_TRFfinds.fasta -evalue 1E-1 -outfmt 7 -db 34KnobsonNCBI_renamed.txt -task blastn -out DB_vincecheck




>gi|31879374|dbj|AB099947.1| Triticum aestivum DNA, centromeric repeat, clone:CS14
CTGGCCTTGAGAAGACGTTCGAAACAAAGCTCGATAACAGATTTAATGAACTGCTTGCGCGTCTTTCACC
ACCGGCTGCACCTGCCGCACCTCTGCAGCAGCAGCAACAACAACAACAACAACAACAACAACAACAACAA
CAACAACAACAAAGACGACTACCTCCACGTCGCGAAACAGCCATCCGTCGAGCGAGCCGTGTCCTTCTTT
AGCCTGGCCAAACTGTTGGTGCTGCTGTTAATACTTTTGTGGCTCCTGCTGCTGATGCGGAGGATGATTA
TGCGGGAGATTACGAGG


>gi|31879374|dbj|AB099947.1| Triticum aestivum DNA, centromeric repeat, clone:CS14
CTGGCCTTGAGAAGACGTTCGAAACAAAGCTCGATAACAGATTTAATGAACTGCTTGCGCGTCTTTCACCACCGGCTGCACCTGCCGCACCTCTGCAGCAGCAGCAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAAAGACGACTACCTCCACGTCGCGAAACAGCCATCCGTCGAGCGAGCCGTGTCCTTCTTTAGCCTGGCCAAACTGTTGGTGCTGCTGTTAATACTTTTGTGGCTCCTGCTGCTGATGCGGAGGATGATTATGCGGGAGATTACGAGG













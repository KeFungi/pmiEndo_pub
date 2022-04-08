# phylogenetics
## BUSCO
```bash
#wd=busco/
conda install -c bioconda -c conda-forge busco=4.1.1
spcodel=$(cat ../meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
    busco -m transcriptome -i ../CDS/${spcode}_CDS.fasta -o ${spcode}_fungi -l fungi_odb10 --offline
done

```

## align each busco
```bash
#wd=busco/
Rscript ../scripts/filter_busco_gene.R
../scripts/make_busco_fasta.sh busco_fungi_gene_tb.tsv picked_fasta

#wd=busco/picked_fasta
for f in *.fasta; do
  outpath=${f%.fasta}_ali.fasta
  mafft --auto $f > $outpath
done

#wd=busco/
../scripts/concat_busco_fasta.py picked_fasta ../meta_tables/rel_genome_list.txt #write to supermatrix_ali.fasta
```

## RAxML tree from CIPRESS
run raxml with supermatrix_ali.fasta in blackbox mode at http://cipress.com
Saved to ./phylotree/RAxML_busco_65.newick and ./phylotree/RAxML_busco_65_bootstrap.newick


# gene number stat

## genome completeness
```bash
#wd=/
Rscript scripts/genome_completeness.R meta_tables/rel_genome_list.txt meta_tables/genome_completeness.tsv
```
## genome size
```bash
#wd=/
Rscript scripts/genome_size.R meta_tables/rel_genome_list.txt meta_tables/genome_size.tsv
```
## gene count
```bash
#wd=/
Rscript scripts/genome_gene_count.R meta_tables/rel_genome_list.txt meta_tables/genome_gene_count.tsv
```
## gene length
```bash
#wd=/
Rscript scripts/genome_gene_length.R meta_tables/rel_genome_list.txt meta_tables/genome_gene_length.tsv
```

# SignalP
```bash
#wd=signalp
spcodel=$(cat ../meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
  signalp -fasta ../proteins/${spcode}.aa.fasta -org euk -format short -prefix $spcode
done
```

# TargetP
```bash
#wd=targetp
spcodel=$(cat ../meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
  targetp -fasta ../proteins/${spcode}.aa.fasta -org non-pl -format short -prefix $spcode
done
```

# WoLFPSort
```bash
#wd=WoLFPSort
spcodel=$(cat ../meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
  runWolfPsortSummary fungi < ../proteins/${spcode}.aa.fasta > ${spcode}.WoLFPSort
done
```

# tmhmm
```bash
#wd=tmhmm
spcodel=$(cat ../meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
  tmhmm --short ../proteins/${spcode}.aa.fasta > ${spcode}.tmhmm
done
```

# ps_scan
```bash
#wd=ps_scan
wget ftp://ftp.expasy.org/databases/prosite/prosite.dat
spcodel=$(cat ../meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
  perl ps_scan.pl -d prosite.dat -p PS00014 ../proteins/${spcode}.aa.fasta > ${spcode}_PS00014.scan
done
```

# secretome pipe
```bash
#move everything to /secretome
cp targetp secretome; cp  secretome; cp WoLFPSort secretome; cp signalp secretome; cp tmhmm secretome

#wd=/
Rscript scripts/secretome_pipe.R #write summary to secretome/secretome_count.tsv
```

# parse repeat
```bash
#wd=/
mkdir repeat
rm results/repeat_length.tsv

spcodel=$(cat meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
  zcat ge/${spcode}.fasta.gz | perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (sort(keys(%fa))){ print "$s\n$fa{$s}\n" }}' | perl -lne 'if(/^>(\S+)/){ $n=$1} else { while(/([a-z]+)/g){ printf("%s\t%d\t%d\n",$n,pos($_)-length($1),pos($_)) } }' > results/${spcode}.mask.bed
  sum=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' repeat/${spcode}.mask.bed)
  echo ${spcode}$'\t'${sum} >> results/repeat_length.tsv
done
```

# CAZymes
https://mycocosm.jgi.doe.gov/mycocosm/annotations/browser/cazy/summary?p=fungi
```bash
#wd=cazyme
wget 'https://mycocosm.jgi.doe.gov/mycocosm/annotations/browser/cazy/summary?p=fungi' -O Cazyme_table.html
../scripts/html2csv.py Cazyme_table.html Cazyme_table.csv
```

# SMC: Secondary metabolite clusters (also where genome names are from)
https://mycocosm.jgi.doe.gov/pages/sm-clusters-summary.jsf?organism=fungi&genomes=fungi
copy the full table as .csv by hand as *JGI_SMC_raw.csv*
use MS excel to extract USR from strain name hyperlink (which include portal id)
```bash
#wd=SMC
wget 'https://mycocosm.jgi.doe.gov/pages/sm-clusters-summary.jsf?organism=fungi&genomes=fungi' -O JGI_SMC_raw.csv
Rscript parse_JGI_SMC.R #extract portal id from URL; write to JGI_SMC.csv
```

# pfam
From mycocosm https://mycocosm.jgi.doe.gov/mycocosm/annotations/browser/pfam/summary?p=fungi
```bash
#wd=pfam
taxa=$(cat ../meta_tables/mycocosm_main_taxa.txt)
mkdir raw_pfam

for i in $taxa
do
  wget "https://mycocosm.jgi.doe.gov/mycocosm/annotations/browser/pfam/summary?p=${i}" -O raw_pfam/${i}_pfam.html #download pfam page
  grep '<table' raw_pfam/${i}_pfam.html > raw_pfam/${i}_pfam.html.grep #extract table
  ../scripts/html2csv.py raw_pfam/${i}_pfam.html.grep raw_pfam/${i}_pfam.csv
done

#wd=raw_pfam
Rscript ../concat_csv.R #write to /pfam/pfam.csv
```

# RGI
```bash
#wd=rgi
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json ./card.json --local

spcodel=$(cat ../meta_tables/rel_genome_list.txt)

for spcode in $spcodel
do
    rgi main --input_sequence ../ge/${spcode}.fasta --output_file $spcode --input_type contig --clean --local
done

#wd=/
Rscript scripts/parse_RGI.R
```

# protease
```bash
merops_db=diamond/merops.dmnd
mkdir merops
rm results/merops.tsv
spcodel=$(cat meta_tables/rel_genome_list.txt)
for spcode in $spcodel
do
  diamond blastp -p 1 -d $merops_db -q proteins/${spcode}.aa.fasta -o merops/${spcode}.merops.tsv
  sum=$(awk '{print $1}' merops/${spcode}.merops.tsv | sort | uniq | wc -l)
  echo ${spcode}$'\t'${sum} >> results/merops.tsv
done
```

# lipase
```bash
led_db=diamondled.dmnd
mkdir led
rm results/led.tsv
spcodel=$(cat meta_tables/rel_genome_list.txt)
for spcode in $spcodel
do
  diamond blastp -p 1 -d $led_db -q proteins/${spcode}.aa.fasta -o led/${spcode}.led.tsv
  sum=$(awk '{print $1}' led/${spcode}.led.tsv | sort | uniq | wc -l)
  echo ${spcode}$'\t'${sum} >> results/led.tsv
done
```

# statistical tests and plots

```bash
Rscript statistics_and_plots.R
```

# CAZyme GH88 DTL
1. all GH88 from DTL/GH88_count.csv downloaded as GH88.fasta
2. global alignment mafft-eisi

exclude short sequences
```
jgi|Sodal1|333603|fgenesh1_kg.6_#_98_#_Locus6055v2rpkm4.77
jgi|LecAK0013_1|411950|fgenesh1_kg.15_#_42_#_TRINITY_DN8590_c0_g1_i1
jgi|Conlig1|706301|MIX18391_10_51
jgi|Botci1|7364|BC1T_05486
jgi|Cadsp1|635582|MIX34793_2_35
```
save to DTL/GH88_curate.fasta

3. update DTL/GH88_count_curate.cvs
   mafft-einsi GH88_curate.fasta > GH88_curate_ali.fasta
   convert GH88_curate_ali.fasta to GH88_curate_ali.fasta.phy

4. tree build by RAxML at CIPRESS
```
MLsearch_CAT_ false
convergence_criterion_  false
datatype_ protein
disable_ratehet_  false
disable_seqcheck_ false
intermediate_treefiles_ false
mesquite_output_  false
mulcustom_aa_matrices_  false
no_bfgs_  false
outsuffix_  result
parsimony_seed_val_ 12345
printbrlength_  false
prot_matrix_spec_ DAYHOFF
prot_sub_model_ PROTCAT
provide_parsimony_seed_ true
rearrangement_yes_  false
runtime_  0.25
select_analysis_  fd
specify_nchar_  1000
specify_runs_ false
```
output DTL/GH88.tree

5. midpoint root to GH88.newick
6. parse by DTL/ranger_input.R
  ouput DTL/rager_input.newick
7. run Ranger-DTL.mac -i rager_input.newick -o output.reconciliation
   The minimum reconciliation cost is: 170 (Duplications: 3, Transfers: 46, Losses: 26)



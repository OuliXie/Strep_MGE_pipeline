# MGE scripts

Scripts used to extract and classify accessory genomic segments in Streptococcus dysgalactiae subsp. equisimilis (SDSE)/S. pyogenes transmission study in remote communities in the NT, Australia (https://doi.org/10.1101/2023.08.17.23294027).
Automated MGE class calling from pangenomes using Corekaburra (https://github.com/milnus/Corekaburra) and Panaroo (https://github.com/gtonkinhill/panaroo) outputs based on recombinase types within accessory segments.
This is an update of the pipeline described at https://github.com/OuliXie/Global_SDSE

Modified from the MGE classification scheme developed by Khedkar S. et al. NAR 2022 DOI: 10.1093/nar/gkac163.
Requires recombinase profile HMMs available at https://promge.embl.de/download.cgi, Smyshlyaev G. et al. Mol Sys Biol 2021 (https://doi.org/10.15252/msb.20209880), and pfam (https://www.ebi.ac.uk/interpro/) with pfam profile HMM names provided in Supplementary Table 1 in the ProMGE publication by Khedkar et al.
Requires type 4 secretion system (T4SS) profile HMMs from MacSyFinder (https://github.com/gem-pasteur/Macsyfinder_models/tree/master/models/TXSS/profiles) doi:10.1371/journal.pone.0110726

A convenience script is provided to download the profile HMMs

These scripts are tested with the following dependencies and versions
- Panaroo v1.2.10
- Corekaburra v0.0.5
- hmmer 3.3.2
- Biopython v1.79
- Python v3.7.12
- Numpy v1.21.6
- Pandas v1.3.4
- csvtk v0.23.0
- samtools v1.15.1
- bedtools v2.30.0
- pybedtools v0.9.0
- emapper v2.1.7 (with eggNOG DB v5.2.0)

## Updates from previous version described at https://github.com/OuliXie/Global_SDSE

- Improved accessory segment extraction algorithm by saving contig position in ordered gene_presence_absence files (prevents extracting across contigs)
- Extracts fasta and gff of all accessory segments from each genome
- Post-processing script provided to summarise MGEs to user-specified insertion sites

## Usage

Should first construct pangenome using Panaroo and analyse pangenome gene synteny using Corekaburra. Gffs used in Corekaburra and this pipeline should be updated after running Panaroo to correct for filtered and refound genes (https://gtonkinhill.github.io/panaroo/#/post/output_gffs).
Classification of accessory segments will be influenced by the diversity of genomes used to build the pangenome. e.g., if all isolates are from a closely related outbreak, a universally present prophage will be classified as core and will be missed by this pipeline.
Detection of IS may also be influenced by parameters used to construct pangenome. For example, Panaroo in `--clean-mode strict` will tend to remove IS annotations at the end of contigs and may result in underestimation of IS elements (particularly for SDSE). 
Should run Panaroo with `-a pan` to output alignment of all genes.

**Step 1.**
Download profile HMMs.
Run `$bash get_hmms.sh`

Will download required profile HMMs from proMGE, pfam and MacSyFinder.
Profile HMMs will be placed in subdirectories `hmms/recombinase` and `hmms/T4SS`

Need to manually download "Source Data for Appendix" from Smyshlyaev G. et al. Mol Sys Biol 2021 (https://doi.org/10.15252/msb.20209880), extract zip file and move `TRdb.hmm` to `hmms/recombinase`

**Step 2.**
Run `$ python get_pangenome_protein.py -i <path/to/directory/containing/MSA>`

Translates nucleotide sequence from Panaroo pangenome to protein sequence and gives a representative pangenome.
Takes a folder of multisequence alignments and/or multifasta containing sequences of a single gene e.g., MSA from Panaroo.
Excludes sequences with an internal stop codon and takes longest remaining sequence per gene.
Generates representative protein pangenome file pan_genome_reference_protein.fa

**Step 3.**
Run `$ bash mge_hmmsearch.sh -g <pangenome_reference_protein.fa> -i <gene_presence_absence.csv> -p <path/to/hmms/folder>`

Runs hmmsearch against pangenome genes to annotate recombinase/integrase genes from mobile genetic elements (MGEs).
Looks for recombinase subfamilies, phage structural proteins and essential type IV secretion system proteins (coupling protein, ATPase).
Uses eggnog-mapper v2.1.7 DIAMOND search under --sensitive mode for annotating phage structural proteins.
Consider running eggnog-mapper separately as takes >30 mins - commented out of script by default

e.g. `$ python emapper.py -i pan_genome_reference_protein.fa -o sdse`

Expects hmms folder to contain subdirectories hmms/recombinase and hmms/T4SS containing their profile hmms

**Step 4.**
Run `$ python order_gene_presence_absence.py -i <annotated_gene_presence_absence.csv> -g <path/to/postpanaroo_gffs/directory> -o <output> (default ordered_gene_presence_absence)`

Generates gene_presence_absence_roary file ordered by sequence ID for each sequence.
Uses the order of CDS from post-Panaroo corrected gffs to order the genes for each genome.
This is because sorting by locus tag can result in refound genes being out of order.
Required for pulling out accessory genes for each core gene pair.
Expects annotated_gene_presence_absence.csv file which has been annotated by mge_hmmsearch.sh.

**Step 5.**
`filter_core_pair_summary.sh` is the main script used to generate analysis outputs.
```
$ bash filter_core_pair_summary.sh 
  -n <max_acc lower threshold>
  -i <path/to/corekaburra/output/folder>
  -r <recombinase_rules.tsv>
  -g <path/to/gffs/folder>
  -f <path/to/all/genome/assemblies_fasta/folder>
  -c <core threshold> (default 0.99)
  -m <min accessory genes per genome if pass -n cutoff> (default 1)
  -a <min accessory genes in coreless fragments> (default 10)
  -b <ignore gene classifications from sequence breaks>
```

Requires as input a folder containing all assemblies for genomes used in the pangenome `-f` and a folder containing all gffs `-g`. These files are required to allow automated extraction of accessory segment fasta and gffs. All accessory segments are extracted regardless of MGE classification - this allows analysis of non-MGE genes in addition to MGEs.

Note that `-n` are `-m` are different. 
`-n` selects core-core pairs to investigate which have at least one genome with `-n` accessory genes in the segment (taken from max_acc in Corekaburra).
`-m` determines which genomes to include when investigating a core-core pair based on an accessory gene cutoff between the core-core pair.
For example, `-n 1` will select core-core pairs which have a max_acc of >=1. 
`-m 1` will analyse segments which have >=1 accessory genes between core-core pairs which passed `-n`.
This combination of parameters may be used to focus analysis on accessory segments of sufficient size to carry cargo genes which may be of biological interest but will also allow detection of insertion sequences which have only one CDS at these insertion sites. Suggest using `-n 1` and `-m 1` for most comprehensive analysis.

Beware, setting `-n` to a value >1 may cause a small number of genes to be classified incorrectly as MGE-related when they are actually non-MGE e.g., when non-MGE gene lies of an accessory segment containing only 1 gene but a rearrangement might exist when a MGE inserts next to it.

By default, will only look at accessory only segments/contigs (coreless) with >= 10 accessory genes.
This is to balance between finding missing elements and not including fragmented partial MGEs/poor assemblies.
Have found that smaller thresholds e.g., 5 accessory genes, when looking at accessory only segments/contigs may lead to calling of poorly supported small elements.
Can change this cutoff using `-a` flag.
If interested in small plasmids, may need to drop this threshold. Note that plasmids will be binned with "ICE" elements - these scripts were developed for analyses of S. pyogenes and S. dysgalactiae subsp. equisimilis (SDSE) which uncommonly carry plasmids. Could search within "coreless" segments for elements classified as ICE which may be plasmids in other organisms.

`-b` enhances accessory gene classification and is recommended. When using rules to determine class of gene (e.g., non-MGE, phage etc.), occurrences of the gene in accessory segments with a sequence break which are not classified as a MGE, will not be weighted. These frequently occur in the setting of a sequence break where accessory genes are on the other side of a sequence break from MGE recombinases. These occurrences will instead be binned in an "Unclassified category".

Calls:
- per_pair_accessory.py
- coreless_extract.py
- mge_type.py
- summarise_mge_count.py
- summarise_gene_type.py

Generates 5 summary outputs:
- *summary_acc_count.tsv* <br/>
a tab delimited file listing number of accessory genes in a segment with genome labels as rows and core-core pairs as columns
- *summary_dis_acc_summary.tsv* <br/>
similar to summary_acc_count.tsv except with accessory segment length instead of gene count
- *summary_mge_count.csv* <br/>
a comma delimited file listing classification of MGE type (or non-MGE) for each accessory segment including element count for nested MGEs in the same format as summary_acc_count (e.g., nested element with 3 IS recombinases and 1 phage will be labelled IS&3:Phage&1).
Also includes a "Hotspot" designation when >=4 recombinases are present in a single element or when both phage and ICE structural genes are present in an accessory segment to flag possible nested elements.
- *summary_mge_classes.csv* <br/>
similar to *summary_mge_count.csv* but with element counts removed
- *classified_genes.csv* <br/>
a comma delimited file in the style of Roary gene_presence_absence.csv with the addition of classification of each genes into core, non-MGE accessory, or MGE type. In addition, the frequency each gene appears within a MGE type is listed. 

Extracts accessory segment fasta and gffs:
- fastas are saved in a subdirectory designated by their respective flanking core gene pair within `per_pair_output/segment_fasta/` <br\>
fastas are extracted from the nucleotide after the first flanking core gene to the nucleotide immediately before the second flanking core gene or the end of a contig in the case of a sequence break
- gffs are saved in a subdirectory designated by their respective flanking core gene pair within `per_pair_output/segment_gff/` <br\>
gffs are extracted from the first accessory CDS to the last accessory CDS

Additionally, interim files with the accessory segment and MGE classifications for each genome and core-core combination are saved in a subdirectory `per_pair_output/` for manual inspection if desired.
Accessory only segments (coreless) are saved to dummy location tags e.g., coreless_1, coreless_2 etc.
A list of which coreless contigs they belong to is saved as `coreless_segments.csv`.

The `per_pair_output/` directory contains a large number of files. Consider saving the `segment_fasta`, `segment_gff` and `coreless_contigs` subdirectories but otherwise removing the rest of this directory to free up space post-analysis and after inspection.

# Post-processing scripts

```
$ python order_mge.py
  -i <csv of user defined insertion regions>
  -m <summary_mge_classes.csv>
  -f <path to segment_fasta folder with all accessory segment fastas>
  -o <output> (default ordered_segments>
```

Extracts all segment fastas categorised as phage, phage_like, ICE or ME (including if nested with IS).
Takes as input a two-column comma-separated file with user-defined insertion regions (using pangenome names) which can be given a new alias. This file should NOT have headings. The first column designates the pair of core genes forming an insertion region separated by a hyphen `-`. The second column is an optional alias for this insertion region. If it is left blank, the two core genes designating the insertion region will be used.

For example:

| eno-sagA | 20 |
| ahpC-Sequence_break | 52 |

If given two core genes (e.g., eno-sagA), the script will automatically summarse any segment with a contig break (Sequence_break) but adjacent to one of these core genes (e.g., eno or sagA) to the insertion region.
In the above example, any MGE between eno-sagA, eno-Sequence_break, Sequence_break-eno, sagA-Sequence_break, Sequence_break-sagA will be grouped with eno-sagA.
These will be output to a subdirectory within the specified output directory named by the alias in the second column of the user given csv (e.g., in the case of eno-sagA, it will be called i_20). 
Therefore, only the most common arrangement between the core genes should be used else genome re-arrangements with contig breaks may be grouped with this insertion region incorrectly. 
However, if the most common arrangement between core genes is used and frequency of rearrangements are checked (this can be confirmed using Corekaburra's `core_pair_summary.tsv` output), incorrect classification should be uncommon. 
Fully assembled rearrangements will NOT be summarised (e.g., eno-rpsA is a rearrangement and will not be grouped with eno-sagA).

Any insertion regions containing phage, phage_like, ICE or ME not specified in the user given csv will still be extracted but will not be summarised as above and will be output into a subdirectory given by its flanking core gene pair (or Sequence_break).

Fastas will be named \<genome\>_\<insertion\>.fa and fasta headers will contain the genome name, contig location and MGE classification delimited by `;`.

# Additional information
**per_pair_accessory.py, mge_type.py, and coreless_extract.py**
```
$ python per_pair_accessory.py
  -f <core gene 1>
  -s <core gene 2>
  -l <low_frequency_gene_placement.tsv>
  -r <recombinase_rules.tsv>
  -o <output directory (default per_pair_output/)>
  -fa <directory of fasta files>
  -gf <directory of gff files>
  -i <sequence_id.txt>
  -g <path/to/ordered_gene_presence_absence/directory>
  -c <core threshold> (default 0.99)
  -m <min accessory genes per genome if pass -n cutoff> (default 1)
```
```
$ python mge_type.py
  -i <accessory segment csv extracted by per_pair_accessory.py>
  -r <recombinase_rules.tsv>
  -o <output directory>
  -b <ignore gene classification if no MGE found - use only when sequence break in gene pair>
```
```
$ python coreless_extract.py
  -a <Corekaburra's coreless_contig_accessory_gene_content.tsv>
  -f <path/to/gffs/directory>
  -fa <path/to/fastas/directory>
  -o <output directory (default per_pair_output/)>
  -i <sequence_id.txt>
  -g <path/to/ordered_gene_presence_absence/directory>
  -c <core threshold> (default 0.99)
  -m <min num,ber of genes in accessory fragment cutoff> (default 10) 
 ```

per_pair_accessory.py, coreless_extract.py and mge_type.py are called by filter_core_pair_summary.sh.

per_pair_accessory.py extracts accessory segments between core gene 1 and core gene 2.
The order of these need to be specified as the script initially assumes gene 1 occurs before gene 2 in the gene_presence_absence.csv.
Also extracts fasta and gff of accessory segments (regardless of if they're classified as MGE or not).

coreless_extract.py finds accessory segments with >= 10 accessory genes. 
A threshold of 10 is used to balance between finding missing elements and not including poorly assembled small segments.
Have found that smaller thresholds e.g., 5, when looking at accessory only segments/contigs may lead to calling of poorly supported small elements.
Will all extract fasta and gff of accessory segments (regardless of if they're classified as MGE or not).

mge_type.py will look for presence of annotated recombinase, T4SS and recombinase genes and classify if the accessory segment belongs to a MGE and which type.

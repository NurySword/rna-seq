All scripts can be found in this folder:
https://github.com/NurySword/rna-seq/blob/main/Scripts/ALL%20SCRIPTS.txt

Where did you get your normal/controls from? Same PRJN or a separate PRJN?
Normal and controls are from the same NCBI project: PRJNA751582. 

What is the URL of your data? What was the original format? How did you convert it to fastqc acceptable format?
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA751582&o=acc_s%3Aa
SRA Tool fasterq-dump has been used which downloaded data in fastq format from NCBI Sequence Read Archive:

fasterq-dump --split-files SRR15326622 -p
fasterq-dump --split-files SRR15326623 -p
fasterq-dump --split-files SRR15326624 -p
fasterq-dump --split-files SRR15326625 -p
fasterq-dump --split-files SRR15326626 -p
fasterq-dump --split-files SRR15326627 -p

Did she/he clean up all her sequences using fastqc/flexbar/fastp? (5 points)
FastQC tool has been utilized to check the data quality of all six SRA data sets hence six data quality check report files has been provided under QC folder:
SRR15326622_fastqc.html
SRR15326623_fastqc.html
SRR15326624_fastqc.html
SRR15326625_fastqc.html
SRR15326626_fastqc.html
SRR15326627_fastqc.html
Also MultiQC tool has been used to generate and view consolidated report for all six data sets. MultiQS report file is also provided under QC folder:
Multiqc_report.html

Flexbar has been used to remove adapter trims.

Did she/he index/align her genome sequences successfully with HISAT2?

331958 reads; of these:
  331958 (100.00%) were paired; of these:
    81564 (24.57%) aligned concordantly 0 times
    244498 (73.65%) aligned concordantly exactly 1 time
    5896 (1.78%) aligned concordantly >1 times
    ----
    81564 pairs aligned concordantly 0 times; of these:
      22742 (27.88%) aligned discordantly 1 time
    ----
    58822 pairs aligned 0 times concordantly or discordantly; of these:
      117644 mates make up the pairs; of these:
        70545 (59.96%) aligned 0 times
        45113 (38.35%) aligned exactly 1 time
        1986 (1.69%) aligned >1 times
89.37% overall alignment rate


331958 reads; of these:
  331958 (100.00%) were paired; of these:
    79998 (24.10%) aligned concordantly 0 times
    246080 (74.13%) aligned concordantly exactly 1 time
    5880 (1.77%) aligned concordantly >1 times
    ----
    79998 pairs aligned concordantly 0 times; of these:
      23327 (29.16%) aligned discordantly 1 time
    ----
    56671 pairs aligned 0 times concordantly or discordantly; of these:
      113342 mates make up the pairs; of these:
        69625 (61.43%) aligned 0 times
        41532 (36.64%) aligned exactly 1 time
        2185 (1.93%) aligned >1 times
89.51% overall alignment rate


331956 reads; of these:
  331956 (100.00%) were paired; of these:
    80827 (24.35%) aligned concordantly 0 times
    245108 (73.84%) aligned concordantly exactly 1 time
    6021 (1.81%) aligned concordantly >1 times
    ----
    80827 pairs aligned concordantly 0 times; of these:
      23447 (29.01%) aligned discordantly 1 time
    ----
    57380 pairs aligned 0 times concordantly or discordantly; of these:
      114760 mates make up the pairs; of these:
        71059 (61.92%) aligned 0 times
        41283 (35.97%) aligned exactly 1 time
        2418 (2.11%) aligned >1 times
89.30% overall alignment rate


390607 reads; of these:
  390607 (100.00%) were paired; of these:
    96183 (24.62%) aligned concordantly 0 times
    287099 (73.50%) aligned concordantly exactly 1 time
    7325 (1.88%) aligned concordantly >1 times
    ----
    96183 pairs aligned concordantly 0 times; of these:
      30454 (31.66%) aligned discordantly 1 time
    ----
    65729 pairs aligned 0 times concordantly or discordantly; of these:
      131458 mates make up the pairs; of these:
        85236 (64.84%) aligned 0 times
        44368 (33.75%) aligned exactly 1 time
        1854 (1.41%) aligned >1 times
89.09% overall alignment rate


390607 reads; of these:
  390607 (100.00%) were paired; of these:
    96029 (24.58%) aligned concordantly 0 times
    287368 (73.57%) aligned concordantly exactly 1 time
    7210 (1.85%) aligned concordantly >1 times
    ----
    96029 pairs aligned concordantly 0 times; of these:
      31392 (32.69%) aligned discordantly 1 time
    ----
    64637 pairs aligned 0 times concordantly or discordantly; of these:
      129274 mates make up the pairs; of these:
        86234 (66.71%) aligned 0 times
        41121 (31.81%) aligned exactly 1 time
        1919 (1.48%) aligned >1 times
88.96% overall alignment rate


390607 reads; of these:
  390607 (100.00%) were paired; of these:
    96366 (24.67%) aligned concordantly 0 times
    286957 (73.46%) aligned concordantly exactly 1 time
    7284 (1.86%) aligned concordantly >1 times
    ----
    96366 pairs aligned concordantly 0 times; of these:
      30711 (31.87%) aligned discordantly 1 time
    ----
    65655 pairs aligned 0 times concordantly or discordantly; of these:
      131310 mates make up the pairs; of these:
        88862 (67.67%) aligned 0 times
        40373 (30.75%) aligned exactly 1 time
        2075 (1.58%) aligned >1 times
88.63% overall alignment rate


BAM STATISTICS FOR NORMAL (NORMAL.BAM)

2056990 + 0 in total (QC-passed reads + QC-failed reads)
65246 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
1845761 + 0 mapped (89.73% : N/A)
1991744 + 0 paired in sequencing
995872 + 0 read1
995872 + 0 read2
1506966 + 0 properly paired (75.66% : N/A)
1688288 + 0 with itself and mate mapped
92227 + 0 singletons (4.63% : N/A)
37816 + 0 with mate mapped to a different chr
34769 + 0 with mate mapped to a different chr (mapQ>=5)


BAM STATISTICS FOR TREATED (TUMOR.BAM)

2419653 + 0 in total (QC-passed reads + QC-failed reads)
76011 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
2159321 + 0 mapped (89.24% : N/A)
2343642 + 0 paired in sequencing
1171821 + 0 read1
1171821 + 0 read2
1766486 + 0 properly paired (75.37% : N/A)
1995542 + 0 with itself and mate mapped
87768 + 0 singletons (3.74% : N/A)
13282 + 0 with mate mapped to a different chr
11297 + 0 with mate mapped to a different chr (mapQ>=5)


Did she/he finish StringTie runs and created alignment bams successfully?
Yes. Please see below.
Did she/he run StringTie to create the FPKM counts?
Yes. All head and tail section outputs of normal and treated (tumor) files can be found in this google drive folder (coverage-fpkm-tpm). 



Did she/he run StringTie to create the raw counts? 
Yes. All head and tail section outputs of normal and treated (tumor) files can be found in this google drive folder.


Did she/he run Ballgown to identify DE genes and create visualizations (including heatmap)?
All R related script files and Rout files can be found in this folder (R scripts and Rout).
PDF output files can be found here.

A sample plot is below:


Did she/he run rnaseqmut to create the vcf files?
Yes.

Did she/he run annovar to annotate the vcf results and pick the important variations? 

Did she/he run a gene set enrichment analysis on both DE genes and annovar genes? 


Gene ID
Gene Definition / Impact
2687
This gene is a member of the gamma-glutamyl transpeptidase gene family, and some reports indicate that it is capable of cleaving the gamma-glutamyl moiety of glutathione. The protein encoded by this gene is synthesized as a single, catalytically-inactive polypeptide, that is processed post-transcriptionally to form a heavy and light subunit, with the catalytic activity contained within the small subunit. The encoded enzyme is able to convert leukotriene C4 to leukotriene D4, but appears to have distinct substrate specificity compared to gamma-glutamyl transpeptidase. Alternative splicing results in multiple transcript variants encoding different isoforms. [provided by RefSeq, Oct 2014]
Expression
Ubiquitous expression in fat (RPKM 18.6), kidney (RPKM 15.2) and 24 other tissues
9127
The protein encoded by this gene belongs to the family of P2X receptors, which are ATP-gated ion channels and mediate rapid and selective permeability to cations. This gene is predominantly expressed in skeletal muscle, and regulated by p53. The encoded protein is associated with VE-cadherin at the adherens junctions of human umbilical vein endothelial cells. Alternative splicing results in multiple transcript variants. A related pseudogene, which is also located on chromosome 22, has been identified. [provided by RefSeq, Apr 2009]
706
Present mainly in the mitochondrial compartment of peripheral tissues, the protein encoded by this gene interacts with some benzodiazepines and has different affinities than its endogenous counterpart. The protein is a key factor in the flow of cholesterol into mitochondria to permit the initiation of steroid hormone synthesis. Alternatively spliced transcript variants have been reported; one of the variants lacks an internal exon and is considered non-coding, and the other variants encode the same protein. [provided by RefSeq, Feb 2012]
2687
This gene is a member of the gamma-glutamyl transpeptidase gene family, and some reports indicate that it is capable of cleaving the gamma-glutamyl moiety of glutathione. The protein encoded by this gene is synthesized as a single, catalytically-inactive polypeptide, that is processed post-transcriptionally to form a heavy and light subunit, with the catalytic activity contained within the small subunit. The encoded enzyme is able to convert leukotriene C4 to leukotriene D4, but appears to have distinct substrate specificity compared to gamma-glutamyl transpeptidase. Alternative splicing results in multiple transcript variants encoding different isoforms. [provided by RefSeq, Oct 2014]
29775
The caspase recruitment domain (CARD) is a protein module that consists of 6 or 7 antiparallel alpha helices. It participates in apoptosis signaling through highly specific protein-protein homophilic interactions. Like several other CARD proteins, CARD10 belongs to the membrane-associated guanylate kinase (MAGUK) family and activates NF-kappa-B (NFKB; see MIM 164011) through BCL10 (MIM 603517) (Wang et al., 2001 [PubMed 11259443]).[supplied by OMIM, Mar 2008]
85358
This gene is a member of the Shank gene family. Shank proteins are multidomain scaffold proteins of the postsynaptic density that connect neurotransmitter receptors, ion channels, and other membrane proteins to the actin cytoskeleton and G-protein-coupled signaling pathways. Shank proteins also play a role in synapse formation and dendritic spine maturation. Mutations in this gene are a cause of autism spectrum disorder (ASD), which is characterized by impairments in social interaction and communication, and restricted behavioral patterns and interests. Mutations in this gene also cause schizophrenia type 15, and are a major causative factor in the neurological symptoms of 22q13.3 deletion syndrome, which is also known as Phelan-McDermid syndrome. Additional isoforms have been described for this gene but they have not yet been experimentally verified. [provided by RefSeq, Mar 2012]
23654
Members of the B class of plexins, such as PLXNB2 are transmembrane receptors that participate in axon guidance and cell migration in response to semaphorins (Perrot et al. (2002) [PubMed 12183458]).[supplied by OMIM, Mar 2008]
7078
This gene belongs to the TIMP gene family. The proteins encoded by this gene family are inhibitors of the matrix metalloproteinases, a group of peptidases involved in degradation of the extracellular matrix (ECM). Expression of this gene is induced in response to mitogenic stimulation and this netrin domain-containing protein is localized to the ECM
706
Present mainly in the mitochondrial compartment of peripheral tissues, the protein encoded by this gene interacts with some benzodiazepines and has different affinities than its endogenous counterpart. The protein is a key factor in the flow of cholesterol into mitochondria to permit the initiation of steroid hormone synthesis. Alternatively spliced transcript variants have been reported; one of the variants lacks an internal exon and is considered non-coding, and the other variants encode the same protein. [provided by RefSeq, Feb 2012]



Did she/he investigate her/his Gene Set Enrichment Analysis (GSEA) results and explain during the demo (15 points)

Did she/he write the report to explain the findings? (10 points)
Below genes were my focus genes based on expression difference analysis performed by SupplementaryR.R file execution. 

Changes in gene expressions suggest:
changes in proinflammatory T-cell marker system (genes: CARD10, IGLV3-10, NCF4, APOBEC3D)
Changes in glucose metabolism (genes: GGT1, ACO2, 

Gene ID
Gene Definition/Impact
CARD10

increase
CARD10 (Caspase Recruitment Domain Family Member 10) is a Protein Coding gene. Diseases associated with CARD10 include Immunodeficiency 89 And Autoimmunity and Glaucoma, Primary Open Angle. Among its related pathways are TCR Signaling (Qiagen) and GPCR Pathway. Gene Ontology (GO) annotations related to this gene include signaling receptor complex adaptor activity. An important paralog of this gene is CARD11.
IGLV3-10

decrease
V region of the variable domain of immunoglobulin light chains that participates in the antigen recognition (PubMed:24600447). Immunoglobulins, also known as antibodies, are membrane-bound or secreted glycoproteins produced by B lymphocytes. In the recognition phase of humoral immunity, the membrane-bound immunoglobulins serve as receptors which, upon binding of a specific antigen, trigger the clonal expansion and differentiation of B lymphocytes into immunoglobulins-secreting plasma cells. Secreted immunoglobulins mediate the effector phase of humoral immunity, which results in the elimination of bound antigens (PubMed:20176268, PubMed:22158414). The antigen binding site is formed by the variable domain of one heavy chain, together with that of its associated light chain. Thus, each immunoglobulin has two antigen binding sites with remarkable affinity for a particular antigen. The variable domains are assembled by a process called V-(D)-J rearrangement and can then be subjected to somatic hypermutations which, after exposure to antigen and selection, allow affinity maturation for a particular antigen (PubMed:17576170, PubMed:20176268). 
LV310_HUMAN,A0A075B6K4
NCF4

decrease
The protein encoded by this gene is a cytosolic regulatory component of the superoxide-producing phagocyte NADPH-oxidase, a multicomponent enzyme system important for host defense. This protein is preferentially expressed in cells of myeloid lineage. It interacts primarily with neutrophil cytosolic factor 2 (NCF2/p67-phox) to form a complex with neutrophil cytosolic factor 1 (NCF1/p47-phox), which further interacts with the small G protein RAC1 and translocates to the membrane upon cell stimulation. This complex then activates flavocytochrome b, the membrane-integrated catalytic core of the enzyme system. The PX domain of this protein can bind phospholipid products of the PI(3) kinase, which suggests its role in PI(3) kinase-mediated signaling events. The phosphorylation of this protein was found to negatively regulate the enzyme activity. Alternatively spliced transcript variants encoding distinct isoforms have been observed. [provided by RefSeq, Jul 2008]


APOBEC3D

decrease
This gene is a member of the cytidine deaminase gene family. It is one of seven related genes or pseudogenes found in a cluster, thought to result from gene duplication, on chromosome 22. Members of the cluster encode proteins that are structurally and functionally related to the C to U RNA-editing cytidine deaminase APOBEC1. The protein encoded by this gene catalyzes site-specific deamination of both RNA and single-stranded DNA. The encoded protein has been found to be a specific inhibitor of human immunodeficiency virus-1 (HIV-1) infectivity. [provided by RefSeq, Mar 2017]


GGT1

decrease
The enzyme encoded by this gene is a type I gamma-glutamyltransferase that catalyzes the transfer of the glutamyl moiety of glutathione to a variety of amino acids and dipeptide acceptors. The enzyme is composed of a heavy chain and a light chain, which are derived from a single precursor protein. It is expressed in tissues involved in absorption and secretion and may contribute to the etiology of diabetes and other metabolic disorders. Multiple alternatively spliced variants have been identified. There are a number of related genes present on chromosomes 20 and 22, and putative pseudogenes for this gene on chromosomes 2, 13, and 22. [provided by RefSeq, Jan 2014]


GAS2L1

increase
This gene encodes a member of the growth arrest-specific 2 protein family. This protein binds components of the cytoskeleton and may be involved in mediating interactions between microtubules and microfilaments. This protein localizes to the proximal end of mature centrioles and links centrosomes to both microtubules and actin. Alternate splicing results in multiple transcript variants. A pseudogene of this gene is found on chromosome 9. [provided by RefSeq, May 2018]


PRAME

increase
This gene encodes an antigen that is preferentially expressed in human melanomas and that is recognized by cytolytic T lymphocytes. It is not expressed in normal tissues, except testis. The encoded protein acts as a repressor of retinoic acid receptors, and likely confers a growth advantage to cancer cells via this function. Alternative splicing results in multiple transcript variants. [provided by RefSeq, Apr 2014]


ACO2

decrease
ACO2 (Aconitase 2) is a Protein Coding gene. Diseases associated with ACO2 include Infantile Cerebellar-Retinal Degeneration and Optic Atrophy 9. Among its related pathways are Pyruvate metabolism and Citric Acid (TCA) cycle and Glucose metabolism. Gene Ontology (GO) annotations related to this gene include iron ion binding and aconitate hydratase activity. An important paralog of this gene is ACO1.
EMID1

decrease
Predicted to be located in several cellular components, including Golgi apparatus; endoplasmic reticulum; and extracellular matrix. Predicted to be part of collagen trimer. [provided by Alliance of Genome Resources, Nov 2021]



Currently, no oral medications are available for individuals suffering from type 1 diabetes (T1D). Our randomized placebo-controlled phase 2 trial recently revealed that oral verapamil has short- term beneficial effects in subjects with new-onset type 1 diabetes (T1D) 1. However, what exact biological changes verapamil elicits in humans with T1D, how long they may last, and how to best monitor any associated therapeutic success has remained elusive. We therefore now conducted extended analyses of the effects of continuous verapamil use over a 2-year period, performed unbiased proteomics analysis of serum samples and assessed changes in proinflammatory T-cell markers in subjects receiving verapamil or just standard insulin therapy. In addition, we determined the verapamil-induced changes in human islets using RNA sequencing. Our present results reveal that verapamil regulates the thioredoxin system and promotes an anti-oxidative and anti-apoptotic gene expression profile in human islets, reverses T1D-induced elevations in circulating proinflammatory T-follicular-helper cells and interleukin-21 and normalizes serum levels of chromogranin A (CHGA), a recently identified T1D autoantigen 2,3. In fact, proteomics identified CHGA as the top serum protein altered by verapamil and as a potential therapeutic marker. Moreover, continuous use of oral verapamil delayed T1D progression, promoted endogenous beta cell function and lowered insulin requirements and serum CHGA levels for at least 2 years and these benefits were lost upon discontinuation. Thus, the current studies provide crucial mechanistic and clinical insight into the beneficial effects of verapamil in T1D. Overall design: Treatment of islet cells with compound (verapamil) followed by transcriptomics. 

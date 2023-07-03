# 2. Sequence Data

**Preamble: What Data Are We Using?**

!!! info ""

    **Genome In A Bottle & HG002** 

    In this workshop we will be using data from HG002, which is a reference sample from the [Genome In A Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) consortium. The GIAB project releases benchmark data for genomic characterization, and you may have seen their benchmark variant calls and regions out in the wild. As part of their benchmarking material generation, they release datasets for their reference samples. We will be using those in this workshop.
    
    **Family Structure**
    
    HG002 is actually part of a trio of reference samples. Below is the family listing, also known as the Ashkenazim trio:
    
    * HG002: Son
    * HG003: Father
    * HG004: Mother
    
    If you'd like to use data from any of the Ashkenazim trio, there is a set of [index files on the GIAB github](https://github.com/genome-in-a-bottle/giab_data_indexes).
    
    **There is an excellent HG002 assembly available**
    
    Since HG002 has so much data to compare against, it is pretty common for new technologies to benchmark on HG002&mdash;making even more HG002 data. It has grown into a data ecosystem. This is part of the reason why it was chosen for the T2T consortium as the target for one of its next big pushes: a [high-quality diploid T2T genome](https://github.com/marbl/HG002). This assembly has been worked on by many leading people in the field. And while it is undergoing polishing and is talked about as a draft, it is gapless outside of rDNA arrays and is very good. This makes HG002 a very good sample for testing assembly processes. You can always go back and compare your results against the draft from the T2T consortium.
    
## Data For Graph Building: PacBio & ONT Data

We are going to start by introducing the two long read sequencing technologies that we will be using: PacBio's HiFi and Oxford Nanopore's Ultralong. These two technologies are complementary and each have their own strengths. You can answer some questions more easily with HiFi and some more easily with ONT UL. They can also be used together, and this is important for the concept of a hybrid assembly algorithm where accurate reads are used to create a draft assembly and long reads are used to extend that assembly. 

In this section, we will learn about both technologies and then create plots showing their characteristic read lengths and qualities. This will help us get a feel for what the data actually looks like in the wild. Lastly, we will prepare this data for use in our assembly section. As someone who is about to make an assembly, you have the most control over what type of data you put into the assembly algorithm. The more you know about the data, the better your assembly will be.


## PacBio Hifi: Illumina-Like Quality With Long Reads
**What is PacBio HiFi**<br>
PacBio's high fidelity (or HiFi) reads are long (~15kb) and accurate (~99.9%). PacBio produces such high quality reads (with their single-molecule real-time, or SMRT, sequencing) by reading the same sequence over and over again in order to create a circular consensus sequence (or CCS) as shown below. 

**PacBio's CCS Process**
![PacBio CCS Process](https://raw.githubusercontent.com/human-pangenomics/hprc-tutorials/GA-workshop/assembly/genomics_aotearoa/images/sequencing/HiFi-reads-img.svg)

Long, highly accurate reads allows for a number of analyses that were difficult or impossible in the context of short reads. For instance, variants can be more easily phased as you can just look for variants that are seen in the same sequencing read since HiFi read lengths span much more than short read lengths. In our context, long accurate reads allow assembly algorithms to build assembly graphs across difficult regions. But it turns out that HiFi reads aren't long enough to span exact repeats in regions like human centromeres.

## ONT Ultralong: Lower Quality But Really Long Reads
Oxford Nanopore's ultralong (UL) sequencing has lower accuracy (~97%), but is really long (even longer than normal ONT). This is achieved though a different library prep&mdash;as compared to normal DNA sequencing with ONT. UL library prep uses a transposase to cut DNA at non-specific sites where it can then be adapted for sequencing. 

**ONT Ultralong Library Prep**
<p align="center">
    <img src="https://raw.githubusercontent.com/human-pangenomics/hprc-tutorials/GA-workshop/assembly/genomics_aotearoa/images/sequencing/ULK114_workflow_V1-3.svg" width="350"/>
</p>

The time, transposase amount, and temperature are all factors that affect transposase activity. The more you cut, the shorter the reads. ONT's standard DNA library prep, on the other hand, shears DNA then ligates adapters. (If you've created DNA libraries using Illumina's TruSeq kits, then you get the idea.)

These UL reads, while less accurate than HiFi, span tricky regions, which makes UL and HiFi data highly complementary, especially in the context of *de novo* assembly, as we will soon see.


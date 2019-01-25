# GOWSH

Perl homology searcher based on webscrapping and heuristic approaches. It's supposed to look up in HomoloGene and
Ensemble Compara after running Bidirectional best hit algorithm (BDBH) and [OrthoMCL algorithm](https://www.ncbi.nlm.nih.gov/pubmed?Db=pubmed&Cmd=ShowDetailView&TermToSearch=12952885) (OMCL).

## Getting Started

Clone the repo on local:

    git clone https://github.com/carrascomj/gowsh

## Usage

gowsh.pl is the main script. The program takes command-line arguments with
the following options.

    gowsh.pl --gfile|go|glist "path_to_file|GOid|list" --tfile|torg "path_to_file|organism"
        [--modelf|modelo] "path_to_file|organism"

        --gfile path_to_file: input, genes as multiFASTA
        --go GOid: input, Genetic Ontology ID (as in AmiGO)
        --glist list: input, blank separated list gene IDs
        --tfile path_to_file: multiFASTA containg proteins of genome of target organism
        --torg organism: target organism name (genus and specie)
        --modfile path_to_file: optional, multiFASTA containg proteins of genome of model organism
        --modorg organism: optional, model organism name (genus and specie)

## Running the test

The script can be tested from an input file on 'examples' directory.

    gowsh.pl --gfile examples/input_genes.faa --torg 'corynebacterium glutamicum' --modorg 'corynebacterium efficiens'

The program will then parse the input file, download both genomes from NCBI and try to match homologues.

## What I Learned

This code was developed as a project for one subject of my BSc in Biotechnology (UPM). To sum up, I learnt the following concepts:
* Webscrapping GO using Perl.
* Use of Entrez [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/) programmatic access API from NCBI.
* Use of [Ensembl](http://www.ensembl.org/index.html) API.
* Run BLAST on local using [blast+](https://www.ncbi.nlm.nih.gov/pubmed/20003500?dopt=Citation).
* Heuristic algorithms to account for homology.
* How to write a README.md.

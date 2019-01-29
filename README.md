<p align="center">
  <img src="logo.svg" alt="logo" width="400"/>
</p>

# GOWSH

Perl homology searcher based on webscrapping and heuristic approaches. It's supposed to look up in HomoloGene and
Ensemble Compara after running Bidirectional best hit algorithm (BDBH) and [Inparanoid](http://software.sbc.su.se/cgi-bin/request.cgi?project=inparanoid) (OMCL).

## Getting Started

Clone the repo on local:

    git clone https://github.com/carrascomj/gowsh

Add script to path (on your bash initialization file; e.g., ~/bashrc):

    export PATH=$PATH:"path/to/gowsh/bin"

The program requires additional packages that can be installed with [cpanm](https://metacpan.org/pod/cpanm), if not already done:

    cpanm JSON Data::Dumper Bio::SeqIO LWP::Simple File::Basename Getopt::Long XML::Parser

Alternatively, one could install WebAPIsGOWSH as an usual perl package (on 'gowsh/' directory):

    perl Makefile.PL
    make
    make install

Finally, blastall and blast+ are both required.

## Usage

gowsh.pl is the main script. The program takes command-line arguments with
the following options:

    gowsh.pl --gfile|go|glist "path_to_file|GOid|list" --tfile|torg "path_to_file|organism"
        [--modelf|modelo] "path_to_file|organism" --out "outfile" --preserve

        --gfile path_to_file: input, genes as multiFASTA
        --go GOid: input, Genetic Ontology ID (as in AmiGO)
        --glist list: input, blank separated list gene IDs
        --tfile path_to_file: multiFASTA containing proteins of genome of target organism
        --torg organism: target organism name (genus and specie)
        --modfile path_to_file: optional, multiFASTA containing proteins of genome of model organism
        --modorg organism: optional, model organism name (genus and specie)
        --out "outfile": optional, name of output file; default "GOWSH_output.txt"
        --preserve: optional, if it's added, (nearly) all files generated will be preserved.

## Running the test

The script can be tested from an input file on 't/' directory.

    gowsh.pl --gfile t/input_genes.faa --torg 'corynebacterium glutamicum' --modorg 'corynebacterium efficiens'

Otherwise, you can

The program will then parse the input file, download both genomes from NCBI and try to match homologues.

## What I Learned

This code was developed as a project for one subjects of my BSc in Biotechnology (UPM). To sum up, I learnt the following concepts:
* Webscrapping biological information using Perl and [mygene API](http://mygene.info/v3/api#/).
* Use of Entrez [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/) programmatic access API from NCBI.
* Use of [Ensembl REST](http://www.ensembl.org/index.html) API.
* Run BLAST on local using [blast+](https://www.ncbi.nlm.nih.gov/pubmed/20003500?dopt=Citation).
* Heuristic algorithms to account for homology.
* How to buil a Perl package.
* How to write a README.md.

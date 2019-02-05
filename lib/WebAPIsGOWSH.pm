#!/usr/bin/perl

# Package that contains Web APIs to be used by main script

package WebAPIsGOWSH;

use warnings;
use strict;
use Bio::SeqIO;
use LWP::Simple;
use Data::Dumper;
use JSON;

use vars qw(@ISA @EXPORT);

BEGIN {
   require Exporter;
   $WebAPIsGOWSH::VERSION = '0.1';
   @ISA = qw(Exporter);
   @EXPORT = qw(@path_tmp @path_files $out_file %out_table mygeneAPI d_genes d_go d_entrez esearch elink efetch search_homologues EnsemblREST prot2homologene get_stream add_output inparanoidWeb)
}

# Global variables
our $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
our @path_tmp; # array of paths of temporary files
our @path_files; # [0] gfile, [1] tfile y [2] modfile
our $out_file; # output file
our %out_table; # results
our $plant; # just to know which Ensemble to use
our %common_names = (            # mygene API accepts 9 common names.
    "homo sapiens" => "human",   # Arabidopsis isn't included here because plants
    "mus musculus" => "mouse",   # receive different threatment.
    "rattus novergicus" => "rat",
    "drosophila melanogaster" => "fruitfly",
    "caenorhabditis elegans" => "nematode",
    "danio rerio" => "zebrafish",
    "xenopus tropicalis" => "frog",
    "sus scrofa" => "pig"
);

# 1. WebAPI functions
sub mygeneAPI{
    # API of mygene.info
    # INPUTS -> id: string, ID or gene name,
    #               or Entrez protein -> v.g., 'refseq:NM_001798',
    #               or Ensembl protein -> v.g., 'ensembl.protein:ENSP00000243067',
    #               or ontology -> v.g., 'go:0000307';
    #           spec: string, optional, name of specie to filter the search
    # OUTPUT -> hash (reference), parsed json.
    our %common_names; # memoization with genus and specie
    my $id = $_[0];
    my $spec = $_[1] ? $_[1] : "";
    my $in = $_[1];
    if ($spec){
        if (exists $common_names{$spec}){
            # one of the 9 common names was provided
            $spec = "&specie=".$common_names{$spec}
        } else{
            # it will be looked up on NCBI and stored on common_names
            $spec =~ s/ /+/;
            $spec = esearch($spec, "taxonomy");
            $common_names{$in} = $spec;
            $spec = "&specie=".$spec; # memoization
        }
    }
    my @ids_out;
    my %hit;
    if (! $id){
        warn "Llamada a mygeneAPI sin argumentos"
    } elsif ($id =~ /^\d+$/ || $id =~ /^ENS/){
        # Entrez or Ensembl ID
        %hit = get("http://mygene.info/v3/gene/$id" . "?fields=ensembl.gene,entrezgene&dotfield=True")
    } else{
        # Name of gene
        # top 30 (mygene score based)
        my $res_json = get("http://mygene.info/v3/query?q=$id"."&fields=ensembl.protein,refseq.protein,ensembl.gene,symbol,homologene.genes&size=30&dotfield=True$spec");
        %hit = %{ decode_json($res_json) };
    }
    if (%hit){
        return \%hit
    } else{
        return 0
    }
}

sub d_genes{
    # Function to look for genes at mygene. It returns Entrez and Ensembl IDs.
    # INPUTS -> all_genes: array of genes/proteins.
    #           modorg: string, model organism (to filter).
    my @all_genes = @{ $_[0] };
    my $modorg = $_[1];
    my @ids;
    foreach my $gene (@all_genes) {
        my @hits = &mygeneAPI($gene,$modorg)->{hits};
        my $uid; my $ensid;
        if (ref $hits[0] eq 'ARRAY'){ # JSON decoder fails to parse the output sometimes
            @hits = @{ $hits[0] }
        }
        foreach my $hit (@hits) { #sorted by score
            my %hit_hash = %{ $hit };
            if ($uid && $ensid){
                last # we're just taking the best scored template, maybe it should be changed.
            }
            if (! $uid){
                $uid = exists $hit_hash{"refseq.protein"} ? $hit_hash{"refseq.protein"} : 0;
                $uid = ref $uid eq 'ARRAY'? @$uid[0] : $uid;
            } if (! $ensid){
                $ensid = exists $hit_hash{"ensembl.protein"} ? $hit_hash{"ensembl.protein"} : 0;
                $ensid = ref $ensid eq 'ARRAY'? @$ensid[0] : $ensid;
                if ($ensid =~ /mrna$/i){
                    $ensid = 0
                }
            }
        }
        push @ids, $uid, $ensid
    }
    return @ids
}

sub d_go{
    # Function that returns genes from a Gene Ontology.
    # INPUT -> GO id.
    # OUTPUT -> Entrez and Ensembl proteins IDs.
    my $go = "go:".$_[0];
    my @ids;
    my %res = %{ &mygeneAPI($go) };
    my @hits = @{ $res{hits} };
    foreach my $hit (@hits) { # sorted by score
        my %hit_hash = %{ $hit };
        foreach my $entry (qw/refseq.protein ensembl.protein/) {
            my $id = exists $hit_hash{$entry} ? $hit_hash{$entry} : 0;
            if (ref $id eq 'ARRAY'){
                push @ids, @{ $id }
            } elsif ($id){
                push @ids, $id
            }
        }
    }
    return @ids
}
sub d_entrez{
    # Function that downloads a genome or proteins of a genome from NCBI to multiFASTA
    # making use of Eutils API. It respects the 3 calls/s restriction.
    # - input -> query_name: string, organism to find.
    #           dbto: string; usually "aa", "nucleotide", "homologene"
    #           dbfrom: string, where to look up the UID (usually "gene" o "genome")
    #           orgmod: model organism
    # - output -> string, path to downloaded file.
    our @path_tmp;
    our @path_files;
    my $query_name = $_[0];
    my $dbto = $_[1];
    my $dbfrom = $_[2] ? $_[2] : "genome";
    my $orgmod = $_[3];
    print "Searching NCBI-$dbfrom...\n";
    if ($orgmod){
        $orgmod =~ s/ /+/;
        $orgmod = '+AND+' . $orgmod . '[organism]'# genus and specie need to be separated by "+"
    } else{
        $orgmod = ''
    }
    $query_name =~ s/ /+/;
    my $suff = ".faa";
    if ($dbto =~ /a|(?:aminoacid)|(?:protein)/i){
        $dbto = "protein";
    } elsif($dbto =~ /(?:nucleotide)|(?:nuccore)/i){
        $dbto = "nuccore";
        $suff = ".fna"
    }
    # 1. ESEARCH. Find UID of genome from query name.
    my $ids = &esearch($query_name, $dbfrom, $orgmod);
    if (! $ids){
        # Query wasn't found
        return
    }
    my $name_out;
    if ($dbto !~ /homologene/){
        # 2. ELINK. Connect UIDs among databases.
        my ($query_key, $webenv) = &elink($dbfrom, $dbto, $ids);
        # 3. EFETCH. Download file.
        $query_name =~ s/\+/_/; # just to beautify output file name
        $name_out = $query_name.$suff;
        &efetch($name_out, $dbto, $query_key, $webenv);
    } else{
        # EFETCH.
        $query_name =~ s/\+/_/;
        $name_out = $query_name . "_homologene" . $suff;
        &efetch($name_out, $dbto, $ids);
    }
    if ($dbfrom =~ /genome/){
        push @path_files, $name_out;
    }
    push @path_tmp, $name_out;
    return $name_out
}

sub esearch{
    # ESEARCH Eutils API. It looks up the query to find UIDs associated.
    # INPUTS -> query_name: string, name to find,
    #           dbfrom: string, where to find UID,
    #           orgmod: name of organism to filter search.
    # OUTPUT -> ids: string, comma separated UIDs.
    my $query_name = $_[0];
    my $dbfrom = $_[1];
    my $orgmod = $_[2] ? $_[2] : "";
    our $base;
    our $plant;
    my $url = $base . "esearch.fcgi?db=$dbfrom&term=$query_name$orgmod";
    my $req_xml = get($url);
    if (!$req_xml){
        return 0
    }
    # Parse req_xml and find UID(s)
    my @UID;
    while ($req_xml =~ /<Id>(\d+?)<\/Id>/gs) {
        push @UID, $1
    }
    my $ids = join(',', @UID);
    if ((lc $dbfrom) eq 'taxonomy' && (! defined $plant)){
        # perform efetch without downloading to fastly see if organisms are plants (to Ensembl)
        $url = $base . "efetch.fcgi?db=$dbfrom&id=$ids";
        my $tax_xml = get($url);
        if ($tax_xml =~ /<ScientificName>.+plantae<\/ScientificName>/){
            $plant = 1;
        } else{
            $plant = 0;
        }
    }
    return $ids
}

sub elink{
    # ELINK Eutils API, links UIDs among databases. Output can be UIDs or history.
    # INPUTS -> dbfrom: string, database where UID comes from (usually "gene" o "genome"),
    #           dbto: string, database where UID will be translated,
    #           ids: string, strictly comma separated UIDS,
    #           $use_h: if not false, it will output a pointer to web history.
    # OUTPUTs -> query_key: string, ID en el historial de Entrez de la búsqueda
    #           webenv: string, entorno web
    our $base;
    my $dbfrom = $_[0];
    my $dbto = $_[1];
    my $ids = $_[2];
    my $use_h = $_[3] ? "" : "&cmd=neighbor_history";
    my $url = $base . "elink.fcgi?dbfrom=$dbfrom&db=$dbto$use_h&id=$ids";
    my $req_xml = get($url);
    if (! $req_xml){
        return 0
    }
    if ($use_h){
        my $query_key = $1 if ($req_xml =~ /<QueryKey>(\d+)<\/QueryKey>/);
        my $webenv = $1 if ($req_xml =~ /<WebEnv>(\S+)<\/WebEnv>/);
        return $query_key, $webenv
    } else{
        my @UID;
        while ($req_xml =~ /<Id>(\d+?)<\/Id>/gs) {
            push @UID, $1
        }
        shift @UID; # 1st UID parsed is from query
        my $ids = join(',', @UID);
    }
}

sub efetch{
    # EFETCH Entrez Eutils API, downloads a file given UIDs or, ideally,a query
    # key and a WebEnv.
    # INPUTS -> name: string, output file name
    #           dbto: string, database where UID is found.
    #           query: string, query key/UIDS
    #           webenv: string, Web env
    # OUPUT -> output file name (just for convenience)
    our $base;
    my $name = $_[0];
    my $dbto = $_[1];
    my $query = $_[2];
    my $webenv = $_[3];
    my $url;
    if ($webenv){
        # if 4 args are passed, it will try do download on buffer.
        $url = $base . "efetch.fcgi?db=$dbto&query_key=$query&WebEnv=$webenv&rettype=fasta&retmode=text"
    } else{
        $url = $base . "efetch.fcgi?db=$dbto&id=$query&rettype=fasta&retmode=text"
    }
    my $out = get($url);
    open (OUT, ">$name");
    print OUT $out;
    close OUT;
    # efetch downloads no more than 10000 sequences, so iterations are needed to
    # download a large genome.
    # TODO: check if  $seqcount is a bottleneck. If not, change to pure Perl.
    my $seqcount = `fgrep -c ">" $name`; # number of sequences in the file
    my $counter = 0;
    my $ret_point = 10000;
    while (($seqcount / 10000) !~ /\./ && $counter<10){ # while num sequences is multiple of 10000
        my $new_url = $url . "&retstart=$ret_point";
        $out = get($new_url);
        open (OUT, ">>$name");
        print OUT $out;
        close OUT;
        $seqcount = `fgrep -c ">" $name`;
        $ret_point += 10000;
        # 100.000 sequences is enough
        $counter++;
    }
    return $name;
}

sub search_homologues{
    # Function that looks up homologues in HomoloGene, Ensembl homlogy and Inparanoid database.
    # INPUTS -> stream_look: Bio::SeqIO stream object,
    #           org_target: string, target organism name,
    our @path_tmp;
    my $stream_run = $_[0];
    my $org_target = $_[1];
    my @prot_ids;
    my $out;
    my $conthg = 0;
    my $contens = 0;
    my $continp = 0;
    seek($stream_run->_fh, 0, 0); # pointer to seq 1
    while(my $seqobj = $stream_run->next_seq()){
        # 1. Refseq IDs
        my $query = $seqobj->id;
        push @prot_ids, $query
    }
    my $not_added_hg_temp = 1;
    print "\nSearching HomoloGene, Ensembl and InParanoid...\n";
    foreach my $prot(@prot_ids){
        # 1. HomoloGene
        my $hgprots_file = &prot2homologene($prot, "temp.faa");
        if ($hgprots_file){ # homologues were found
            if ($not_added_hg_temp){
                push @path_tmp, $hgprots_file;
                $not_added_hg_temp = 0
            }
            my $matches = &get_stream($hgprots_file);
            while (my $seqobj = $matches->next_seq()){
                my $id_match = $seqobj->id;
                # filter by organism
                my $orgfrom = ($seqobj->description =~ /\[([a-z1-9\. ]+)\]$/i)[0];
                if ($orgfrom){ # sometimes organism isn't encoded on descripction
                    $orgfrom = lc $orgfrom;
                    if ($orgfrom eq $org_target){
                        &add_output($prot, $id_match, "HomoloGene");
                        $conthg++;
                    }
                }
            }
        }
        # 3.2. Ensembl homology
        my $res = mygeneAPI("refseq:$prot")->{hits}[0];
        my %hit = %{ $res }; # all ensembl.genes that corresponds with the protein
        if ($hit{"symbol"}){
            my @ensemblegenes;
            if(ref $hit{"symbol"} eq 'ARRAY'){
                @ensemblegenes = @{ $hit{"symbol"} }
            } else{
                @ensemblegenes = ($hit{"symbol"})
            }
            foreach my $ens (@ensemblegenes) {
                # look for homologues
                my $data = &EnsemblREST("homology/id",$hit{"symbol"},$org_target)->{"data"};
                my @homologies = $data->[0]{homologies};
                if ((! @homologies) || ! defined $homologies[0]){
                    next
                }
                @homologies = @{ $homologies[0] }; # [0] -> homologies, [1] -> id
                if (@homologies){
                    my @ensprots;
                    foreach my $match(@homologies) {
                        # add refseq protein name, just the 1st
                        $match = $match->{target}{protein_id};
                        my $ensprot = mygeneAPI("ensembl.protein:$match")->{hits}[0]{"refseq.protein"}; # devolvemos la proteína de refseq
                        if (! $ensprot){
                            $ensprot = $match; # Entrez protein equivalent wasn't found.
                            print "Refseq equivalent wasn't found for in mygene for $ensprot. Ensembl ID will be returned.\n"

                        }
                        if (ref $ensprot eq 'ARRAY'){
                            $ensprot = $ensprot->[0]
                        }
                        push @ensprots, $ensprot;
                    }
                    &add_output($prot, join(',', @ensprots), "Ensembl");
                    $contens++
                }
            }
        }
        # 3.3. Inparanoid online
        my @matches_inp;
        if ($hit{"ensembl.protein"}){
            @matches_inp = inparanoidWeb($hit{"ensembl.protein"}, $org_target);
        } elsif ($hit{"symbol"}){
            @matches_inp = inparanoidWeb($hit{"symbol"}, $org_target);
        } elsif ($hit{"_id"}){
            @matches_inp = inparanoidWeb($hit{"_id"}, $org_target);
        }
        if (@matches_inp){
            &add_output($prot, join(',', @matches_inp), "InparanoidWEB");
            $continp++
        }
    }
    print "$conthg matches where produced by HomoloGene\n";
    print "$contens matches where produced by Ensembl\n";
    print "$continp matches where produced by InParanoid online\n"
}

sub EnsemblREST{
    # Ensembl API data access.
    # INPUTS -> endpoint: string, any of https://rest.ensembl.org/
    #                     v.g., "homology/id"
    #           id: id to look up.
    #           target_specie: string[optional], TaxID or a name
    # OUTPUT -> (reference) json hash ("data" : @array)
    our %common_names; # memoization
    our $plant; # if true, Ensembl genomes will be used.
    my $endpoint = $_[0] ? $_[0] : die "(WebAPIsGOWSH) EnsemblREST API require 2 args, None supplied on line:".__LINE__;
    my $id = $_[1] ? $_[1] : die "(WebAPIsGOWSH) EnsemblREST API require 2 args, 1 ($endpoint) supplied on line:".__LINE__;
    my $target_specie = $_[2];
    my $in = $_[2];
    if (! $target_specie){
        $target_specie = "";
    } elsif ($target_specie =~ /^\d+$/){
        $target_specie = "target_taxon=$target_specie";
    } else{
        if (exists $common_names{$target_specie}){
            if ($common_names{$target_specie} =~ /^\d+$/){ # memoization
                $target_specie = $common_names{$target_specie};
                $target_specie = "target_taxon=$target_specie"
            } else{
                $target_specie = "target_species=$target_specie";
            }
        } else{
            $target_specie =~ s/ /+/;
            $target_specie = esearch($target_specie, "taxonomy");
            $common_names{$in} = $target_specie; # memoization
            $target_specie = "target_taxon=".$target_specie
        }
        $target_specie =~ s/[\+ ]/_/;
    }
    my $base;
    if (! $plant){
        $base = "https://rest.ensembl.org/"
    } else{
        $base = "https://rest.ensemblgenomes.org/"
    }
    my $url = $base . "$endpoint/$id?$target_specie" . ";content-type=application/json";
    my %res_json = %{ decode_json(get($url)) };
    return \%res_json;
}

sub prot2homologene{
    # Function that takes a protein name (refseq) and returns HomoloGene matches.
    # INPUTS -> prot: string, refseq protein,
    #           out: string, opctional, output file name
    # OUTPUT -> string, output file name
    my $prot = $_[0] ? $_[0] : die "prot2homologene require 1 arg, None supplied on line:".__LINE__;
    my $out = $_[1] ? $_[1] : "temphomogene.faa";
    my $uid_prot = &esearch($prot, "protein"); # prot -> uid prot
    my $uid_gene = &elink("protein", "gene", $uid_prot, "no"); # uid prot -> uid gene
    my $uid_hg = &elink("gene", "homologene", $uid_gene, "no"); # uid gene -> uid homologene
    if ($uid_hg){
        $out = &efetch($out, "homologene", $uid_hg);
    }else{
        return 0;
    }
    # Format output to pure multiFASTA.
    open FH,  '<',  $out if -e $out;
    open FHO, '>', ".temp2homogene.faa" or die "Can't write new file: $!";
    my $com = 0;
    while (<FH>){
        if ($com){
            $com = 0;
            next;
        }
        if($_ =~ /^>/){ # consistent ID
            $_ =~ s/gi\|\d+\|ref\|([_a-zA-Z\d\.]+)\|/$1/;
            print FHO $_
        } elsif($_ !~ /\d+\: HomoloGene/){
            print FHO $_
        } else{
            $com = 1
        }
    }
    close FH;
    close FHO;
    system("mv", ".temp2homogene.faa", "$out");
    return $out;
}

sub inparanoidWeb{
    # Webscrapping of InParanoid database.
    # INPUTS -> query: string, protein,
    #           org: string, organism to filter search.
    # OUTPUTS -> string (match) or 0 (unmatched)
    my $query = $_[1];
    my $org = $_[2] ? $_[2] : "all";
    my $url = "http://inparanoid.sbc.su.se/cgi-bin/gene_search.cgi?id=" . $query ."&idtype=all&all_or_selection=all&specieslist=7&scorelimit=0.05&.submit=Submit+Query&.cgifields=specieslist&.cgifields=idtype&.cgifields=all_or_selection";
    my $res_xml = get($url);
    my @res = split "\n", $res_xml;
    my $this_line;
    my $parse_org;
    my $match;
    my $i = 0;
    my @matches;
    # Parse xml requested.
    foreach (@res){
        if ($parse_org){
            if ($_ =~ /(taxonomy.+>$org)/i){
                $parse_org = 0;
                push @matches, $match;
                $i = 0;
            }
            $i++;
            if ($i > 3){
                $parse_org = 0;
                $i = 0;
            }
        } elsif ($this_line){
            $this_line = 0;
            $parse_org = 1;
            if ($_ =~ /([A-Za-z0-9_\.]+)<\/a/){
                $match = "$1"
            }
        } elsif ($_ =~ /<td class=\"proteinid\">/){
            $this_line = 1;
        }
    }
    my @returned_matches;
    # Return just unique IDs.
    foreach my $al (@matches){
        if (@returned_matches){
            if(! grep $_ eq $al, @returned_matches){
                push @returned_matches, $al
            }
        } else{
            push @returned_matches, $al
        }
    }
    return @returned_matches;
}
# 2. I/O functions
sub add_output {
	# Adds a match to hash where results are stored.
	# INPUTS -> key_id: string, gene ID (query)
    #           match: string, gene ID (found)
    #           source: string, method which provided the match.
    our %out_table;
	my $key_id = $_[0];
    my $match = $_[1];
    my $source = $_[2]; # BDBH, HomoloGene, Ensembl...
	if (exists $out_table{$key_id}{$source}){
        push @{$out_table{$key_id}{$source}}, $match
    } else{
        if ($match){
            $out_table{$key_id}{$source} = [$match]
        } else{
            $out_table{$key_id}{$source} = []
        }
    }
}

sub get_stream {
	# Function that pulls out the sequence contained on a multiFASTA
	# INPUTS -> seqs: string, path to file
	#			type: string, format of file
	# - output -> stream object
	my $seqs =$_[0];
	my $type = $_[1] ? $_[1] : "FASTA"; # FASTA as default
	my $inseq = Bio::SeqIO -> new(
        #-string => $seqs,
        -file => $seqs,
		-format => $type, ); # stream is instanced
	return $inseq;
}

"TRUE STATEMENT SO USE RETURNS TRUE";

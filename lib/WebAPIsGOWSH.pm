#!/usr/bin/perl

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
   @EXPORT = qw(@path_tmp @path_files $out_file %out_table mygeneAPI d_genes d_go d_entrez esearch elink efetch search_homologues EnsemblREST prot2homologene extrae_stream add_output inparanoidWeb)
}

our $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
our @path_tmp; # path a archivos temportales que serán borrados
our @path_files; # [0] gfile, [1] tfile y [2] modfile (de haberlo)
our $out_file; # almacena archivo de output
our %out_table; # almacena el output a lo largo del programa
our $plant; # ensembl funciona diferente con plantas
our %common_names = (            # la API de mygene acepta 9 nombres comunes o Taxonomy IDs
    "homo sapiens" => "human",   # aunque se elimina arabidopsis por la implementación
    "mus musculus" => "mouse",
    "rattus novergicus" => "rat",
    "drosophila melanogaster" => "fruitfly",
    "caenorhabditis elegans" => "nematode",
    "danio rerio" => "zebrafish",
    "xenopus tropicalis" => "frog",
    "sus scrofa" => "pig"
);

# 1. WebAPI functions
sub mygeneAPI{
    # Función que devuelve Ensembl o Entrez (NCBI) IDs de un gen y sus proteínas
    # dado uno de estos IDS o el nombre del gen, mediante la API de mygene.info
    # INPUTS -> id: cadena, ID o nombre de gen,
    #           o también una proteína de Entrez -> v.g., 'refseq:NM_001798',
    #           o una proteína de Ensembl -> v.g., 'ensembl.protein:ENSP00000243067',
    #           o una ontología -> v.g., 'go:0000307';
    #           spec: cadena, opcional, para filtrar la búsqueda en base al nombre del gen.
    # OUTPUT -> hash (referencia), contiene el json devuelto por la query.
    our %common_names; # aprovechamos el hash para hacer memoization con género y especie
    my $id = $_[0];
    my $spec = $_[1] ? $_[1] : "";
    my $in = $_[1];
    if ($spec){
        if (exists $common_names{$spec}){
            # se ha proporcionado uno de los 9 nombres comunes (o memo)
            $spec = "&specie=".$common_names{$spec}
        } else{
            # se busca en NCBI-Taxonomy
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
        # tenemos un id de ENTREZ o Ensembl
        %hit = get("http://mygene.info/v3/gene/$id" . "?fields=ensembl.gene,entrezgene&dotfield=True")
    } else{
        # tenemos el nombre de un gen
        # se cogen los 30 de máxima score
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
    # Función que busca genes en un mygene y devuelve los identificadores de Entrez
    # y Ensembl.
    # INPUTS -> all_genes: array, lista de genes/proteínas.
    #           modorg: cadena, organismo modelo.
    my @all_genes = @{ $_[0] };
    my $modorg = $_[1];
    my @ids;
    foreach my $gene (@all_genes) {
        my @hits = &mygeneAPI($gene,$modorg)->{hits};
        my $uid; my $ensid;
        if (ref $hits[0] eq 'ARRAY'){ # fallo en el parseo de decode_json
            @hits = @{ $hits[0] }
        }
        foreach my $hit (@hits) { # ordenados de mayor a menor score
            my %hit_hash = %{ $hit };
            if ($uid && $ensid){
                last
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
    # Función que devuelve una serie de genes de una ontología genética
    # INPUT -> id de ontología genética
    # OUTPUT -> ids de proteínas de Entrez y Ensembl
    my $go = "go:".$_[0];
    my @ids;
    my %res = %{ &mygeneAPI($go) };
    my @hits = @{ $res{hits} };
    foreach my $hit (@hits) { # ordenados de mayor a menor score
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
    # Función que descarga un genoma o proteoma del ncbi en multiFASTA mediante
    # la API de Eutils. Hace 3 llamadas así que en ningún caso será irrespetuosa
    # con ésta (3 llamadas/s como máximo).
    # - input -> query_name: cadena, nombre del organismo/gen a buscar,
    #           dbto: cadena; normalmente "homologene", "aa", "nucleotide",
    #           dbfrom: cadena, de dónde sacamos el UID (normalmente "gene" o "genome")
    #           orgmod: nombre organismo modelo del que buscamos el gen
    # - output -> cadena, path al archivo descargado
    our @path_tmp;
    our @path_files;
    my $query_name = $_[0];
    my $dbto = $_[1];
    my $dbfrom = $_[2] ? $_[2] : "genome";
    # si $dbfrom es "gene", es mejor precisar de qué organismo queremos sacar el gen
    my $orgmod = $_[3];
    print "Buscando en NCBI-$dbfrom...\n";
    if ($orgmod){
        $orgmod =~ s/ /+/;
        $orgmod = '+AND+' . $orgmod . '[organism]'# género va separado de especie por un "+"
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
    # 1. ESEARCH. Encontramos el UID del genoma asociado al nombre dado.
    my $ids = &esearch($query_name, $dbfrom, $orgmod);
    if (! $ids){
        # Query no encontrada
        return
    }
    my $name_out;
    if ($dbto !~ /homologene/){
        # 2. ELINK. Encontramos el registro en la base de datos escogida ligado al UID.
        my ($query_key, $webenv) = &elink($dbfrom, $dbto, $ids);
        # 3. EFETCH. Descargamos la búsqueda a un archivo en local.
        $query_name =~ s/\+/_/; # para que quede mejor el nombre del archivo
        $name_out = $query_name.$suff;
        &efetch($name_out, $dbto, $query_key, $webenv);
    } else{
        # EFETCH.
        $query_name =~ s/\+/_/; # para que quede mejor el nombre del archivo
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
    # ESEARCH API de Entrez, busca el UID asociado al nombre (query_name) aportado.
    # INPUTS -> query_name: cadena, nombre del organismo/gen a buscar,
    #           dbfrom: cadena, de dónde sacamos el UID (normalmente "gene" o "genome"),
    #           orgmod: nombre organismo modelo del que buscamos el gen
    # OUTPUT -> ids: cadena, UIDs separados por comas
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
    # Parseamos el req_xml para buscar el(los) id(s) del genoma.
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
    # ELINK API de Entrez, relaciona UIDs entre bases de datos. Los output se
    # dan para luego hacer un efetch en buffer (en vez de aportar todos los UIDs).
    # INPUTS -> dbfrom: cadena, de dónde sacamos el UID (normalmente "gene" o "genome"),
    #           dbto: cadena, database en la que buscamos el UID,
    #           ids: cadena, UIDS separados estricatamente por coma
    #           $use_h: si se le pasa un valor, devuelve ids en vez de query_key y
    # OUTPUTs -> query_key: cadena, ID en el historial de Entrez de la búsqueda
    #           webenv: cadena, entorno web
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
        shift @UID; # el primer UID es el de la búsqueda
        my $ids = join(',', @UID);
    }
}

sub efetch{
    # EFETCH API de Entrez, descarga un archivo a donde apuntan los UIDs o, preferiblemente,
    # una query key con su WebEnv.
    # INPUTS -> name: string, nombre de archivo de output
    #           dbto: cadena, database en la que buscamos el UID,
    #           query: cadena, query key ó UIDS separados estricatamente por coma
    #           webenv: cadena, entorno web
    # OUPUT -> name, nombre de archivo de output (es 1 argumento, pero es conveniente)
    our $base;
    my $name = $_[0];
    my $dbto = $_[1];
    my $query = $_[2];
    my $webenv = $_[3];
    my $url;
    if ($webenv){
        # si le pasamos 4 argumentos, interpreta que buscamos la query key
        $url = $base . "efetch.fcgi?db=$dbto&query_key=$query&WebEnv=$webenv&rettype=fasta&retmode=text"
    } else{
        $url = $base . "efetch.fcgi?db=$dbto&id=$query&rettype=fasta&retmode=text"
    }
    my $out = get($url);
    open (OUT, ">$name");
    print OUT $out;
    close OUT;
    # efetch descarga como máx 10000 secuencias, hay que iterar
    # TODO: probar si es realmente un cuello de botella $seqcount al hacerlo con
    # Perl estrictamente. Si no, cambiar por una función en Perl.
    my $seqcount = `fgrep -c ">" $name`; # número de secuencias en fichero
    my $counter = 0;
    my $ret_point = 10000;
    while (($seqcount / 10000) !~ /\./ && $counter<10){ # while num secuencias es múltiplo de 10000
        my $new_url = $url . "&retstart=$ret_point";
        $out = get($new_url);
        open (OUT, ">>$name");
        print OUT $out;
        close OUT;
        $seqcount = `fgrep -c ">" $name`;
        $ret_point += 10000;
        # 100.000 secuencias es suficiente para no eternizar el proceso
        $counter++;
    }
    return $name;
}

sub search_homologues{
    # Función que busca en Homologene (NCBI) genes homólogos a los encontrados en el input
    # INPUTS -> stream_look: objeto Bio::SeqIO de stream, input
    #           org_target: cadena, organismo objetivo,
    our @path_tmp;
    my $stream_run = $_[0];
    my $org_target = $_[1];
    my @prot_ids;
    my $out;
    my $conthg = 0;
    my $contens = 0;
    my $continp = 0;
    seek($stream_run->_fh, 0, 0); # puntero a secuencia 1
    while(my $seqobj = $stream_run->next_seq()){
        # 1. Refseq ids de proteínas
        my $query = $seqobj->id;
        push @prot_ids, $query
    }
    my $not_added_hg_temp = 1;
    print "\nBuscando HomoloGene, Ensembl e InParanoid...\n";
    foreach my $prot(@prot_ids){
        # 2. Obtenemos los uids de HomoloGene y Ensemble ids y dereferenciamos
        # 3. Filtramos por organismo objetivo en cada DB
        # 3.1. HomoloGene
        my $hgprots_file = &prot2homologene($prot, "temp.faa");
        if ($hgprots_file){ # se han encontrado homólogos
            if ($not_added_hg_temp){
                push @path_tmp, $hgprots_file;
                $not_added_hg_temp = 0
            }
            my $matches = &extrae_stream($hgprots_file);
            while (my $seqobj = $matches->next_seq()){
                my $id_match = $seqobj->id;
                # filtramos por organismo
                my $orgfrom = ($seqobj->description =~ /\[([a-z1-9\. ]+)\]$/i)[0];
                if ($orgfrom){ # a veces, no está codificado en la descripción el organismo
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
                # buscamos homólogos en Ensembl
                my $data = &EnsemblREST("homology/id",$hit{"symbol"},$org_target)->{"data"};
                my @homologies = $data->[0]{homologies};
                if ((! @homologies) || ! defined $homologies[0]){
                    next
                }
                @homologies = @{ $homologies[0] }; # [0] es homologies, [1] es id
                if (@homologies){
                    my @ensprots;
                    foreach my $match(@homologies) {
                        # añadimos el nombre de la proteína de Entrez, solo el primero
                        $match = $match->{target}{protein_id};
                        my $ensprot = mygeneAPI("ensembl.protein:$match")->{hits}[0]{"refseq.protein"}; # devolvemos la proteína de refseq
                        if (! $ensprot){
                            $ensprot = $match; # no queda otra que devolver el de Ensembl (o Ensembl Plant)
                            print "No se encontró equivalente a refseq en mygene para $ensprot. Se devuelve el identificador de Ensembl.\n"

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
    print "$conthg secuencias han producido un acierto en HomoloGene\n";
    print "$contens secuencias han producido un acierto en Ensembl\n";
    print "$continp secuencias han producido un acierto en InParanoid online\n"
}

sub EnsemblREST{
    # Función API de acceso a Ensemble
    # INPUTS -> endpoint: cadena, cualquiera de https://rest.ensembl.org/
    #                     v.g., "homology/id"
    #           id: id a buscar
    #           target_specie: cadena[opcinal], puede ser un TaxID o un nombre
    # OUTPUT -> (referencia) hash de json (contiene una sola entrada "data" : @array)
    our %common_names; # aprovechamos el hash para hacer memoization con género y especie
    our $plant; # si verdadero, se buscará en Ensembl genomes
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
    # Función que coge un nombre de proteína (refseq) y devuelve las proteínas homólogas
    # a su gen asociado. Dicho así porque tiene que dar un rodeo por el uid del gen.
    # INPUTS -> prot: cadena, proteína refseq,
    #           out: cadena, opcional, nombre de archivo de output
    # OUTPUT -> cadena, archivo donde se han descargado las proteínas homçologas.
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
    # Entrez descarga de Homologene un multiFASTA con comentarios, por alguna razón.
    open FH,  '<',  $out if -e $out;
    open FHO, '>', ".temp2homogene.faa" or die "Can't write new file: $!";
    my $com = 0;
    while (<FH>){
        if ($com){
            $com = 0;
            next;
        }
        if($_ =~ /^>/){ # id consistente
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
    # Función que hace una request a la página de Inparalog para intentar inferir
    # homología online
    # INPUTS -> query: cadena, proteína
    #           org: cadena, organismo a filtrar
    # OUTPUTS -> cadena (match) o 0 (no ha habido match)
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
    foreach (@res){
        if ($parse_org){
            if ($_ =~ /(taxonomy.+>$org)/i){
                #print $+;
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
	# Añade un match al hash donde se almacenan los resultados
	# INPUTS -> key_id: el ID del gen buscado
    #           match: el ID del gen encontrado
    #           source: por qué método se ha encontrado el match
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

sub extrae_stream {
	# Función que extrae el objeto de stream que contiene
	# todas las secuencias de un multifasta.
	# - input -> seqs (cadena): path al archivo,
	#			 formato (cadena)
	# - output -> objeto stream
	my $seqs =$_[0];
	my $type = $_[1] ? $_[1] : "FASTA"; # estamos leyendo de FASTA
	my $inseq = Bio::SeqIO -> new(
        #-string => $seqs,
        -file => $seqs,
		-format => $type, ); #instancio la stream
	return $inseq;
}

"TRUE STATEMENT SO USE RETURNS TRUE";

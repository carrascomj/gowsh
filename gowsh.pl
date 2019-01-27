#!/usr/bin/perl

# Script para buscar genes homólogos de un genoma en otro.
# Pasos en la implementación:
# 1.1) Input a través de geneID de Ensembl o como MULTIFASTA del modelo.
# 2.1) Input Organismo a buscar como ID de Ensembl o como MULTIFASTA.
# 2.2) Input Organismo a buscar por el nombre científico (webscrapping)
# 3.1) Búsqueda por homología.
# 1.2) Input de los genes del modelo a buscar como GO+organismo (webscrapping)
# 3.2) Busqueda por homologen, ProtCluster.
# 3.3) Solo para plantas, buscar intensivamente los no encontradas a partir de un
#       MAS ya definido.

# Autor: Jorge Carrasco Muriel
# Fecha de creación: 22/01/2019

use warnings;
use strict;
use Bio::SeqIO;
use LWP::Simple;
use Getopt::Long;
use Data::Dumper;
use JSON;

# variables globales
our @path_tmp; # path a archivos temportales que serán borrados
our @path_files; # [0] gfile, [1] tfile y [2] modfile (de haberlo)
our $out_file; # almacena archivo de output
our %out_table; # almacena el output a lo largo del programa
our @all_ids_genes; # por si se llama al principio a la API de mygene
our $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
our $usage = "gowsh.pl busca una serie de genes en el organismo proporcionado
----------------------------------------------------------------------------
USO
----------------------------------------------------------------------------
gowsh.pl --gfile|go|glist path_to_file|GOid|lista --tfile|torg path_to_file|organism
    [--modelf|modelo] path_to_file|organism

- --gfile path_to_file: genes de input a buscar por homología como archivo multiFASTA
- --go GOid: ID de ontología genética a buscar por homología
- --glist GOid: lista separada por espacios de ID de genes/proteínas
- --tfile path_to_file: multiFASTA del proteínas del organismo/genes objetivo
- --torg organism: nombre del organismo (género y especie) objectivo
- --modfile path_to_file: opcional, multiFASTA de proteínas del organismo modelo
- --modorg organism: opcional, nombre del organismo (género y especie) modelo

Por defecto, cogerá path_to_file en ambos casos. \n";

# 1. Main y procesamiento de argumentos
sub main{
    my %opts = &proc_args;
    my $org_mod = "";
    if (exists $opts{"mod"}){
        &runnin_gnomes($opts{gin}, $opts{tar}); # heurística BDBH
        $org_mod = lc(($opts{mod}->next_seq()->description =~ /\[([a-z1-9 ]+)\]$/i)[0]); # nombre de organismo
    } else{
        &runnin_gnomes($opts{gin}, $opts{tar});
    }
    my $org_target = lc(($opts{tar}->next_seq()->description =~ /\[([a-z1-9 ]+)\]$/i)[0]); # nombre de organismo
    # TODO: refactorizar la función siguiente para que cada gen sea evaluado por homologene y emsembl simultáneamente
    &search_homologene($opts{gin}, $org_target, $org_mod); # data mining en homologene
    system("rm", "-r", @path_tmp) if $path_tmp[0];
    our %out_table;
    print Dumper(%out_table);
    &write_tsv(qw/Gene BDBH HomoloGene Ensembl/)
}

sub main2check{
    my %opts = &proc_args;
}

sub proc_args{
    # Función de procesado de argumentos.
    # - output -> hash, input => (stream de Bio::SeqIO) {2:3}
    our $usage;
    our $out_file;
    our @path_files;
    sub hprint{ print "$usage"; exit 0 }
    sub clean_ex{ print "No se proporcionó $_[0]\n --help para ver sus opciones \n." }
    my %opts = ();
    if ($#ARGV+1 == 2){
        # Si solo hay 2 argumentos, se asume que se han introducido dos ficheros.
        $opts{"gin"} = &extrae_stream($ARGV[0]);
        $opts{"tar"} = &extrae_stream($ARGV[1]);
        push @path_files, $ARGV[0], $ARGV[1];
        $out_file = "GONG_output.txt";
        return %opts
    }
    my $genesf = ""; # genes a buscar como archivo
    my @genesl = ""; # ids de genes como lista
    my $go = ""; # genes a buscar a sacar de una id de GO
    my $modelo = ""; # nombre org modelo con el que hacer los match recíprocos en la heurística
    my $modelf = ""; # archivo org modelo con el que hacer los match recíprocos en la heurística
    my $targeto = ""; # nombre de organismo a buscar (genoma)
    my $targetf = ""; # archivo de organismo a buscar (genoma)
    my $output = ""; # archivo de output
    GetOptions ("help|?" => \&hprint,
            "gfile:s" => \$genesf,
            "glist:s{,}"   => \@genesl,
            "go:s"  => \$go,
            "modfile:s" => \$modelf,
            "modorg:s" => \$modelo,
            "tfile:s" => \$targetf,
            "torg:s" => \$targeto,
            "out:s" => \$output);
    # Input de genes a buscar
    if ($genesf){
        # por defecto: cogemos el fichero siempre (aunque haya varios inputs)
        $opts{"gin"} = &extrae_stream($genesf);
        push @path_files, $genesf
    } elsif (scalar(@genesl)!=1){ # siempre coge un atributo (está en experimental)
        our @all_ids_genes; # así luego no se tendrá que buscar dos veces
        shift @genesl; # el primer elemento es un espacio, por alguna razón...
        @all_ids_genes = &d_genes(\@genesl, $modelo);
        my $eids;
        foreach (@all_ids_genes) {
            if ($_ !~ /^ENS/){
                $eids .= "$_," # parseamos uids de Entrez
            }
        }
        chop $eids;
        my $uids .= &esearch($eids, "protein");
        $opts{"gin"} = &extrae_stream(&efetch("genes_list.faa", "protein", $uids));
        push @path_tmp, "gene_list.faa";
        push @path_files, "gene_list.faa";
    } elsif ($go){
        # TODO: VAMOS CON EL INPUT COMO GO
        our @all_ids_genes; # así luego no se tendrá que buscar dos veces
        @all_ids_genes = &d_go($go);
        my $eids;
        foreach (@all_ids_genes) {
            if ($_ !~ /^ENS/){
                $eids .= "$_," # parseamos uids de Entrez
            }
        }
        chop $eids;
        my $uids .= &esearch($eids, "protein");
        $opts{"gin"} = &extrae_stream(&efetch("go_$go.faa", "protein", $uids));
        print $opts{"gin"}; exit 0;
        push @path_tmp, "go_$go.faa";
        push @path_files, "go_$go.faa";
        # $opts{"gin"} = &extrae_stream(&efetch("genes_list.faa", "protein", $uids));
    } else{
        clean_ex("genes a buscar.")
    }
    # Input de organismo objetivo
    if ($targetf){
        $opts{"tar"} = &extrae_stream($targetf);
        push @path_files, $targetf
    } elsif ($targeto){
        $opts{"tar"} = &extrae_stream(&d_entrez($targeto, "protein"))
    }
    # Input de modelo
    if ($modelf){
        $opts{"mod"} = &extrae_stream($modelf);
        push @path_files, $modelf
    } elsif ($modelo){
        $opts{"mod"} = &extrae_stream(&d_entrez($modelo, "protein"))
    }
    # Archivos de output
    $out_file = $output ? $output : "GONG_output.txt";
    return %opts
}

# 2. Manejo I/O
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
sub write_tsv{
    # Función escritora del ouput (TSV) a partir de un hash.
    # INPUT -> cols: array, columnas que tendrá el tsv
    my @cols = @_;
    our %out_table;
    our $out_file;
    open FHO, '>', "$out_file" or die "Can't write new file: $!";
    my $header = "";
    foreach (@cols){
        $header .= "$_\t";
    }
    print FHO "$header\n";
    my $lines = "";
    foreach my $gen (keys %out_table) {
        $lines = "$gen\t";
        foreach my $col (@cols) {
            if (exists $out_table{$gen}{$col}){
                foreach my $term(@{%{$out_table{$gen}}{$col}}){
                    $lines .= "$term,"
                }
                chop $lines; # eliminamos la última coma
                $lines .= "\t"
            } else{
                if ($col ne $cols[0]){ # la columna 1 se ignora
                    $lines .= "FALSE\t"
                }
            }
        }
        $lines .= "\n"
    }
    print FHO $lines;
    close FHO;
}

sub abrir_archivo {
  my $ruta_entrada;
  my $archivo = $_[0];
  if (! $archivo){
    print "No se proporcionaron argumentos.\nIntroduzca la ruta de archivo: ";
    $ruta_entrada = <STDIN>;
    print "\n";
    chomp $ruta_entrada;
  } else{
    $ruta_entrada = $archivo;
  }
  if (! open(FILE, $ruta_entrada)){
    print "No se pudo abrir el archivo: $!\n";
    return 0;
  } else{
    open(FILE, $ruta_entrada);
    my @entrada = <FILE>;
    close(FILE);
    return @entrada;
  }
}
sub es_fasta {
  # La función itera sobre el archivo y no permite dos cabeceros seguidos.
  my @fas_wannabe = @_;
  my $ans = 1;
  if ($fas_wannabe[0] !~ /(?!^>\s)^>/g){
    return 0;
  } elsif ($fas_wannabe[@fas_wannabe-1] =~ /^\>/){
    return 0;
  } elsif ($fas_wannabe[1] =~ /^>/){
    return 0;
  }
  for (my $nl = 2; $nl < @fas_wannabe-1; $nl++) {
    if ($fas_wannabe[$nl] =~ /^>/ and $fas_wannabe[$nl+1] =~ /^>/){
      $ans = 0;
      last;
    }
  }
  return $ans;
}

# 3. Webscrapping
sub mygeneAPI{
    # Función que devuelve Ensembl o Entrez (NCBI) IDs de un gen y sus proteínas
    # dado uno de estos IDS o el nombre del gen, mediante la API de mygene.info
    # INPUTS -> id: cadena, ID o nombre de gen,
    #           o también una proteína de Entrez -> v.g., 'q=refseq:NM_001798',
    #           o una proteína de Ensembl -> v.g., 'ensembl.protein:ENSP00000243067',
    #           o una ontología -> v.g., 'go:0000307';
    #           spec: cadena, opcional, para filtrar la búsqueda en base al nombre del gen.
    # OUTPUT -> hash, contiene el json devuelto por la query.
    my $id = $_[0];
    my $spec = $_[1] ? $_[1] : "";
    if ($spec){
        my %common_names = (    # la API de mygene acepta 9 nombres comunes o Taxonomy IDs
            "homo sapiens" => "human",
            "mus musculus" => "mouse",
            "rattus novergicus" => "rat",
            "drosophila melanogaster" => "fruitfly",
            "caenorhabditis elegans" => "nematode",
            "danio rerio" => "zebrafish",
            "arabidopsis thaliana" => "thale-cress",
            "xenopus tropicalis" => "frog",
            "sus scrofa" => "pig"
        );
        if (exists $common_names{$spec}){
            # se ha proporcionado uno de los 9 nombres comunes
            $spec = "&specie=".$common_names{$spec}
        } else{
            # se busca en NCBI-Taxonomy
            $spec =~ s/ /+/;
            $spec = "&specie=".esearch($spec, "taxonomy");
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
        # se coge el de máxima score
        my $res_json = get("http://mygene.info/v3/query?q=$id"."&fields=ensembl.protein,refseq.protein,ensembl.gene&size=30&dotfield=True$spec");
        $res_json = decode_json($res_json);
        %hit = %{ $res_json };
    }
    return %hit;
}

sub d_genes{
    # Función que busca genes en un mygene y devuelve los identificadores de Entrez
    # y Ensembl
    # INPUTS -> all_genes: array, lista de genes/proteínas.
    #           modorg: cadena, organismo modelo.
    my @all_genes = @{ $_[0] };
    my $modorg = $_[1];
    my @ids;
    foreach my $gene (@all_genes) {
        my %hits = &mygeneAPI($gene,$modorg);
        my $uid; my $ensid;
        foreach my $hit (@{ $hits{hits} }) { # ordenados de mayor a menor score
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
    # INPUT ->
    my $go = "go:".$_[0];
    my @ids;
    my %hits = &mygeneAPI($go);
    foreach my $hit (@{ $hits{hits} }) { # ordenados de mayor a menor score
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
    print "Buscando en NCBI...\n";
    our @path_tmp;
    our @path_files;
    my $query_name = $_[0];
    # TODO: check if dbs are in allowed entrez databases
    my $dbto = $_[1];
    my $dbfrom = $_[2] ? $_[2] : "genome";
    # si $dbfrom es "gene", es mejor precisar de qué organismo queremos sacar el gen
    my $orgmod = $_[3];
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
    my $url = $base . "esearch.fcgi?db=$dbfrom&term=$query_name$orgmod";
    my $req_xml = get($url);
    # Parseamos el req_xml para buscar el(los) id(s) del genoma.
    my @UID;
    while ($req_xml =~ /<Id>(\d+?)<\/Id>/gs) {
        push @UID, $1
    }
    my $ids = join(',', @UID);
}

sub elink{
    # ELINK API de Entrez, relaciona UIDs entre bases de datos. Los output se
    # dan para luego hacer un efetch en buffer (en vez de aportar todos los UIDs).
    # INPUTS -> dbfrom: cadena, de dónde sacamos el UID (normalmente "gene" o "genome"),
    #           dbto: cadena, database en la que buscamos el UID,
    #           ids: cadena, UIDS separados estricatamente por coma
    # OUTPUTs -> query_key: cadena, ID en el historial de Entrez de la búsqueda
    #           webenv: cadena, entorno web
    our $base;
    my $dbfrom = $_[0];
    my $dbto = $_[1];
    my $ids = $_[2];
    my $url = $base . "elink.fcgi?dbfrom=$dbfrom&db=$dbto&cmd=neighbor_history&id=$ids";
    my $req_xml = get($url);
    my $query_key = $1 if ($req_xml =~ /<QueryKey>(\d+)<\/QueryKey>/);
    my $webenv = $1 if ($req_xml =~ /<WebEnv>(\S+)<\/WebEnv>/);
    return $query_key, $webenv
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
    return $name;
}

sub search_homologene{
    # Función que busca en Homologene (NCBI) genes homólogos a los encontrados en el input
    # INPUTS -> stream_look: objeto Bio::SeqIO de stream, input
    #           org_target: cadena, organismo objetivo,
    #           org_model: cadena, organismo modelo
    my $stream_run = $_[0];
    my $org_target = $_[1];
    my $org_mod = $_[2];
    my $query;
    my $out;
    my $cont = 0;
    while(my $seqobj = $stream_run->next_seq()){
        $query = $seqobj->id;
        $out = &d_entrez($query, "homologene", "homologene", $org_mod);
        if (! $out){
            next
        } else{
            push @path_tmp, $out
        }
        # Ahora, se parsea el output obtenido eliminando frases informativas
        &clean_homologene($out);
        while (my $seqtarget = &extrae_stream($out)->next_seq()){
            if ($org_target == lc (($seqtarget->description =~ /\[([a-z1-9 ]+)\]$/i)[0])){
                # Tenemos un homólogo
                $cont++;
                my $match = $seqtarget->id;
                &add_output($query,$match, "BDBH");
            }
        }
    }
    print "$cont secuencias han producido un acierto en HomoloGene\n"
}

sub clean_homologene{
    # Función que parsea el output de homologene para que sea procesado correctamente
    # como MULTIFASTA.
    # INPUT -> $out: archivo a parsear
    my $out = $_[0];
    open FH,  '<',  $out if -e $out;
    open FHO, '>', "temp" or die "Can't write new file: $!";
    while (<FH>){
        if ($_ !~ /^\d+\: HomoloGene/){
            print FHO $_;
        }
    }
    close FH;
    close FHO;
    system("mv", "temp", "$out");
}

# 4. Blast
# TODO: igual sería interesante hacer blast remotos si el usuario no tiene blast+.
sub run_blastmas{
	# Función que realiza un blast+ desde system. Los input son algunos de los
    # del propio blast+.
	# - input -> secuencia (next_seq de SeqIO), file db,
    #           opcionales: e-value, num threads.
	# - output -> ruta al resultado de blast.
	our @path_tmp;
	my $prog = $_[0]->alphabet eq "dna" ? "blastn" : "blastp";
	my $query = $_[0]->seq;
	my $id = $_[0]->id;
	my $db_file = $_[1];
	my $e_value = $_[2] ? $_[2] : 0.00001;
    my $num_threads = $_[3] ? $_[3] : (&get_num_cores - 1);
	my @db_title = "$db_file" =~ /^([a-zA-Z0-9\/_]+)\..+/;
	my $path_now = $db_title[0]."dbtemp";
	# Comprueba si existen los archivos blast de la base de datos, si no, los
    # crea (almacena el path en variable global).
	if (! -d $path_now){
		push @path_tmp, $path_now;
		my $molecule = $prog eq "blastn" ? "-p F" : "-p";
		# podría hacerse en /tmp/
		mkdir $path_now;
		system("cp", "$db_file", "$path_now");
		system("formatdb", "-t", "$db_title[0]", "-i", "$path_now/$db_file", "$molecule");
	}
	# 6.2. Llamada al blast.
	open(FH, '>', "$path_now/querytmp",) or die "$!";
	print FH ">$id\n$query";
	close FH;
	#my @tabular_param = qw/qseqid sseqid length ppos/;
	#system("$prog", "-db", "$path_now/$db_file", "-query", "$path_now/querytmp", "-evalue", "$e_value", "-max_target_seqs", "1", "-outfmt", "10 qseqid sseqid length ppos", "-out","$path_now/outmp");
	system("$prog", "-db", "$path_now/$db_file", "-query", "$path_now/querytmp", "-max_target_seqs", "1", "-outfmt", "10 qseqid sseqid length ppos", "-out","$path_now/outmp");
	return "$path_now/outmp";
	# Re: eliminar archivos temporales en mainprint "$long_cut\n"; if ($i>20){last};
}

sub get_num_cores{
    # Función que busca el número de cores.
    # - output -> int, número de núcleos de procesador  o, si falla, 1.
    # archivo de información de linux.
    my $filename = '/proc/cpuinfo';
    my $num_threads;
    if (open(my $fh, '<', $filename)) {
        while (my $row = <$fh>) {
            chomp $row;
            if ($row =~ /^processor\t: (\d+)$/){
                $num_threads = $+
            }
        }
        # en el archivo se empieza a contar desde cero.
        return $num_threads + 1
    } else {
        warn "No se pudo detectar el número de cores, num_threads = 1";
        return 1
    }
}

# 5. Bidirectional best hit algorithm (BDBH)
sub runnin_gnomes {
    # Función que recorre la stream de SeqIO aportada y la compara con un
    # genoma/base de datos en busca de homólogos
	# - input -> objeto stream (Bio::SeqIO -> new) que correremos,
	#			 objeto stream para extraer la seq si hay match
	# - ouput -> diccionario num_aciertos/conteo
    our %out_table; # como es la primera función a ejecutar, hace la estructura del hash.
	my $stream_run = $_[0];
	my $stream_look = $_[1];
    my %seqs_to_look;
    if ($stream_look){
        while(my $seq = $stream_look->next_seq()){
    		# Metemos las secuencias en un hash con id por clave.
    		my $id_look = $seq->id;
    		$seqs_to_look{$id_look} = $seq;
    	}
    }
	my $cont = 0;
	my $matches;
	my $i = 1;
	while(my $seqobj = $stream_run->next_seq()){
        # Probar diferentes e-values.
		my @all_matches = &abrir_archivo(&run_blastmas($seqobj, $path_files[1]));
		if (! @all_matches){
			next;
		}
        # best_match tendrá (query_id, best_match_id, long alineamiento, porcentaje_similaridad)
		my @best_match = split(",", $all_matches[0]);
        my $long_cut;
        if ($path_files[2]){
            # Organismo modelo proporcionado por usuario.
            my $seq_obj_look = $seqs_to_look{$best_match[1]}; # el seqobj cuyo id es el match_id
            my @all_recp_matches = &abrir_archivo(&run_blastmas($seq_obj_look, $path_files[2]));
            my @best_recp_match = split(",", $all_recp_matches[0]);
    		if ($best_recp_match[1] eq $best_match[0]){
    			$long_cut = $seqobj->length >= $seq_obj_look->length ? $seqobj->length : $seq_obj_look->length;
            } else{
                next;
            }
        } else{
            # organismo no proporcionado, la fiabilidad del algoritmo se empobrece.
            $long_cut = $seqobj->length;
        }
		if ($best_match[2] >= 0.66*$long_cut  and $best_match[3] >= 66){
			# Alineamiento al menos del 66% de la mayor de las secuencias y similaridad de, al menos,
			# 66% -> tenemos un ortólogo (probablemente, a lo mejor, más o menos).
			$cont ++;
			my $self_id = $best_match[0];
			my $id_match = $best_match[1];
            &add_output($self_id, $id_match, "BDBH");
	    }
	}
    print "$cont secuencias han producido un acierto mediante BEST RECIPROCAL MATCHES\n";
    # Finalmente, se reinician los punteros de las streams
    foreach (@_) {
        seek($_->_fh, 0, 0);
    }

}

# Línea para hacer el script tanto módulo como ejecutable ("python-like")
&main() if not caller();

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

# variables globales
our @path_tmp; # path a archivos temportales que serán borrados
our @path_files; # [0] gfile, [1] tfile y [2] modfile (de haberlo)
our $out_file; # almacena archivo de output
our %out_table; # almacena el output a lo largo del programa
our $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
our $usage = "gowsh.pl busca una serie de genes en el organismo proporcionado
----------------------------------------------------------------------------
USO
----------------------------------------------------------------------------
gowsh.pl --gfile|go|glist path_to_file|GOid|lista --tfile|torg path_to_file|organism
    [--modelf|modelo] path_to_file|organism

- --gfile path_to_file: genes de input a buscar por homología como archivo multiFASTA
- --go GOid: ID de ontología genética a buscar por homología
- --glist GOid: lista separada por espacios de ID de genes
- --tfile path_to_file: multiFASTA del proteínas del organismo/genes objectivo
- --torg organism: nombre del organismo (género y especie) objectivo
- --modfile path_to_file: opcional, multiFASTA de proteínas del organismo modelo
- --modorg organism: opcional, nombre del organismo (género y especie) modelo

Por defecto, cogerá path_to_file en ambos casos. \n";

# 1. Main y procesamiento de argumentos
sub main_real{
    my %opts = &proc_args;
    my $org_mod = "";
    if (exists $opts{"mod"}){
        &runnin_gnomes($opts{gin}, $opts{tar}, $opts{mod});
        $org_mod = lc(($opts{mod}->next_seq()->description =~ /\[([a-z1-9 ]+)\]$/i)[0]);
    } else{
        &runnin_gnomes($opts{gin}, $opts{tar});
    }
    my $org_target = lc(($opts{tar}->next_seq()->description =~ /\[([a-z1-9 ]+)\]$/i)[0]);
    &d_entrez($opts{gin}, $org_target, $org_mod);
    system("rm", "-r", @path_tmp) if $path_tmp[0];
}

sub main{
    my %opts = &proc_args;
    #my $org_target = lc $opts{tar}->next_seq()->description;
    my $org_target = lc (($opts{tar}->next_seq()->description =~ /\[([a-z1-9 ]+)\]$/i)[0]);
    my $org_mod = lc (($opts{mod}->next_seq()->description =~ /\[([a-z1-9 ]+)\]$/i)[0]);
    print "\n$org_target\n";
    print "\n$org_mod\n";
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
            "glist:s"   => \@genesl,
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
    } elsif (@genesl){
        $opts{"gin"} = &extrae_stream(&d_entrez(@genesl))
    } elsif ($go){
        $opts{"gin"} = &extrae_stream(&d_go(@genesl))
    } else{
        clean_ex("genes a buscar.")
    }
    # Input de modelo
    if ($modelf){
        $opts{"mod"} = &extrae_stream($modelf);
        push @path_files, $modelf
    } elsif ($modelo){
        $opts{"mod"} = &extrae_stream(&d_entrez($modelo, "protein"))
    }
    # Input de organismos objetivo
    if ($targetf){
        $opts{"tar"} = &extrae_stream($targetf);
        push @path_files, $modelf
    } elsif ($targeto){
        $opts{"tar"} = &extrae_stream(&d_entrez($targeto, "protein"))
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
    my $source = $_[2]; # BRM, HomoloGene, Ensembl...
	if (exists $out_table{$key_id}{$source}){
        push @{$out_table{$key_id}{source}}, $match
    }
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
sub d_entrez{
    # Función que descarga un genoma o proteoma del ncbi en multiFASTA mediante
    # la API de Eutils. Hace 3 llamadas así que en ningún caso será irrespetuosa
    # con ésta (3 llamadas/s como máximo).
    # - input -> query_name: cadena, nombre del organismo/gen a buscar,
    #           dbto: cadena; normalmente "homologene", "aa", "nucleotide",
    #           dbfrom: cadena, de dónde sacamos el UID (normalmente "gene" o "genome")
    #  TODO: investigar diferencia entre "nucleotide y "gene"
    #           orgmod: nombre organismo modelo del que buscamos el gen
    # - output -> cadena, path al archivo descargado
    print "Buscando en NCBI-genome...\n";
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
    print "$ids";
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
    my $orgmod = $_[2];
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
}

sub search_homologene{
    # Función que busca en Homologene (NCBI) genes homólogos a los encontrados en el input
    # INPUTS -> stream_look: objeto Bio::SeqIO de stream, input
    #           org_target: cadena, organismo objetivo,
    #           org_model: cadena, organismo modelo
    my $stream_run = $_[0];
    my $org_target =~ $_[1];
    my $org_mod = $_[2];
    my $query;
    my $out;
    while(my $seqobj = $stream_run->next_seq()){
        $query = $seqobj->id;
        $out = &d_entrez($query, "homologene", "homologene", $org_mod);
        # Ahora, se parsea el output obtenido eliminando frases informativas
        open FH,  '<',  $out if -e $out;
		open FHO, '>', "temp" or die "Can't write new file: $!";
		print FHO "{\n";
		while (<FH>){
            if ($_ !~ /^\d+\: HomoloGene/){
                print FHO $_;
            }
		}
		close FH;
		print FHO "}\n";
		close FHO;
		system("mv", "temp", "$out");
        while (my $seqtarget = &extrae_stream($out)->next_seq()){
            # TODO: esto puede entrar en bucle, depende de cómo inteprete Perl.
            if ($org_target == lc (($seqtarget->description =~ /\[([a-z1-9 ]+)\]$/i)[0])){
                # Tenemos un homólogo
                # TODO: escribir el ouput bien de una vez
                print "smth"
            }
        }
        # El organismo se parsea una vez instancio la stream.
    }
}

sub d_go{
    # Función que descarga una serie de genes a multiFASTA a partir de una GO.
    return
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
	my @db_title = "$db_file" =~ /^([\w]+)\..+/;
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
    my $hash_base = (
        "BDBH" => (0),
        "HomoloGene" => (0),
        "Ensembl" => (0)
    );
	my $stream_run = $_[0];
	my $stream_look = $_[2] ? $_[2] : 0;
    my %seqs_to_look;
    if ($stream_look){
        my %seqs_to_look;
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
        if ($stream_look){
            # Organismo modelo proporcionado por usuario.
            my $seq_obj_look = $seqs_to_look{$best_match[1]}; # el seqobj cuyo id es el match_id
            my @all_recp_matches = &abrir_archivo(&run_blastmas($seq_obj_look, $ARGV[0]));
            my @best_recp_match = split(",", $all_recp_matches[0]);
    		if ($best_recp_match[1] eq $best_match[0]){
    			$long_cut = $seqobj->length >= $seq_obj_look->length ? $seqobj->length : $seq_obj_look->length;
            } else{
                next;
            }
        } else{
            # organismo no proporcionado
            $long_cut = $seqobj->length;
        }
		if ($best_match[2] >= 0.66*$long_cut  and $best_match[3] >= 66){
			# Long alineamiento al menos del 66% de la mayor de las secuencias y similaridad de, al menos,
			# 66% -> tenemos un ortólogo (probablemente, a lo mejor, más o menos).
			$cont ++;
			my $self_id = $best_match[0];
			my $id_match = $best_match[1];
            # TODO: El output tiene que ir en un hash de arrays o en un array de arrays
			&add_output("$self_id\t$id_match\n");
	    }
	}
    print "$cont secuencias han producido un acierto mediante BEST RECIPROCAL MATCHES\n"
}

# Línea para hacer el script tanto módulo como ejecutable ("python-like")
&main() if not caller();

#!/usr/bin/perl

# Script that looks up homologues genes based on user input.

# Author: Jorge Carrasco Muriel
# Data of creation: 22/01/2019

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Data::Dumper;
use JSON;
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0) . '/lib';
use WebAPIsGOWSH;

# variables globales
our $usage = "gowsh.pl looks up homologues genes based on user input.
----------------------------------------------------------------------------
USO
----------------------------------------------------------------------------
gowsh.pl --gfile|go|glist path_to_file|GOid|lista --tfile|torg path_to_file|organism
    [--modelf|modelo] path_to_file|organism --out output --preserve

    --gfile path_to_file: input, genes as multiFASTA
    --go GOid: input, Genetic Ontology ID (as in AmiGO)
    --glist list: input, blank separated list gene IDs
    --tfile path_to_file: multiFASTA containing proteins of genome of target organism
    --torg organism: target organism name (genus and specie)
    --modfile path_to_file: optional, multiFASTA containing proteins of genome of model organism
    --modorg organism: optional, model organism name (genus and specie)
    --out outfile: optional, name of output file; default 'GOWSH_output.txt'
    --preserve: optional, if it's added, (nearly) all files generated will be preserved.\n";

# 1. Main y procesamiento de argumentos
sub main{
    my %opts = &proc_args;
    &runnin_gnomes($opts{gin}, $opts{tar}); # heurística BDBH
    my $org_target = lc(($opts{tar}->next_seq()->description =~ /\[([a-z1-9\. ]+)\]$/i)[0]); # nombre de organismo target
    &search_homologues($opts{gin}, $org_target); # data mining
    if (! $opts{preserve}){
        system("rm", "-r", @path_tmp) if $path_tmp[0];
    }
    &write_tsv(qw/Gene BDBH HomoloGene Ensembl InparanoidWEB/)
}

sub proc_args{
    # Función de procesado de argumentos.
    # - output -> hash, input => (stream de Bio::SeqIO) {2:3}
    our $usage;
    our $out_file;
    our @path_files;
    sub hprint{ print "$usage"; exit 0 }
    sub clean_ex{ print $_[0] . " weren't provided\n --help to print usage\n." }
    my %opts = ();
    if ($#ARGV+1 == 2){
        # Si solo hay 2 argumentos, se asume que se han introducido dos ficheros.
        $opts{"gin"} = &get_stream($ARGV[0]);
        $opts{"tar"} = &get_stream($ARGV[1]);
        push @path_files, $ARGV[0], $ARGV[1];
        $out_file = "GOWSH_output.txt";
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
    my $preserve = 0;
    GetOptions ("help|?" => \&hprint,
            "gfile:s" => \$genesf,
            "glist:s{,}"   => \@genesl,
            "go:s"  => \$go,
            "modfile:s" => \$modelf,
            "modorg:s" => \$modelo,
            "tfile:s" => \$targetf,
            "torg:s" => \$targeto,
            "out:s" => \$output,
            "preserve!" => $preserve);
    # Input de genes a buscar
    if ($genesf){
        # por defecto: cogemos el fichero siempre (aunque haya varios inputs)
        $opts{"gin"} = &get_stream($genesf);
        push @path_files, $genesf
    } elsif (scalar(@genesl)!=1){ # siempre coge un atributo (está en experimental)
        shift @genesl; # el primer elemento es un espacio, por alguna razón...
        my @all_ids_prots = &d_genes(\@genesl, $modelo);
        my $eids;
        foreach (@all_ids_prots) {
            if ($_ !~ /^ENS/){
                $eids .= "$_," # parseamos uids de Entrez
            }
        }
        chop $eids;
        my $uids .= &esearch($eids, "protein");
        $opts{"gin"} = &get_stream(&efetch("genes_list.faa", "protein", $uids));
    } elsif ($go){
        my @all_ids_prots = &d_go($go);
        my $eids;
        foreach (@all_ids_prots) {
            if ($_ =~ /^[XNRMP]{2}_/){
                $eids .= "$_," # parseamos uids de Entrez
            }
        }
        chop $eids;
        my $uids .= &esearch($eids, "protein");
        $opts{"gin"} = &get_stream(&efetch("go_$go.faa", "protein", $uids));
        push @path_tmp, "go_$go.faa";
        push @path_files, "go_$go.faa";
        # $opts{"gin"} = &get_stream(&efetch("genes_list.faa", "protein", $uids));
    } else{
        clean_ex("input genes.")
    }
    # Input de organismo objetivo
    if ($targetf){
        $opts{"tar"} = &get_stream($targetf);
        push @path_files, $targetf
    } elsif ($targeto){
        $opts{"tar"} = &get_stream(&d_entrez($targeto, "protein"))
    }
    # Input de modelo
    if ($modelf){
        $opts{"mod"} = &get_stream($modelf);
        push @path_files, $modelf
    } elsif ($modelo){
        $opts{"mod"} = &get_stream(&d_entrez($modelo, "protein"))
    }
    # Archivos de output
    $out_file = $output ? $output : "GOWSH_output.txt";
    $opts{"preserve"} = $preserve;
    return %opts
}

# 2. Manejo I/O
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
        $lines .= "$gen\t";
        foreach my $col (@cols) {
            if (exists $out_table{$gen}{$col}){
                foreach my $term(@{%{$out_table{$gen}}{$col}}){
                    $lines .= "$term,"
                }
                chop $lines; # eliminamos la última coma
                $lines .= "\t"
            } elsif($col ne $cols[0]){
                $lines .= "FALSE\t"
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
    print "Arguments weren't provided.\nWrite path to file: ";
    $ruta_entrada = <STDIN>;
    print "\n";
    chomp $ruta_entrada;
  } else{
    $ruta_entrada = $archivo;
  }
  if (! open(FILE, $ruta_entrada)){
    print "Can't write to file: $!\n";
    return 0;
  } else{
    open(FILE, $ruta_entrada);
    my @entrada = <FILE>;
    close(FILE);
    return @entrada;
  }
}


# 3. Bidirectional best hit algorithm (BDBH)
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
    print "Running BDBH...\n";
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
			next
		}
        # best_match tendrá (query_id, best_match_id, long alineamiento, porcentaje_similaridad)
		my @best_match = split(",", $all_matches[0]);
        my $long_cut;
        if ($path_files[2]){
            # Organismo modelo proporcionado por usuario.
            my $seq_obj_look = $seqs_to_look{$best_match[1]}; # el seqobj cuyo id es el match_id
            my @all_recp_matches = &abrir_archivo(&run_blastmas($seq_obj_look, $path_files[2]));
            if (@all_recp_matches){
                my @best_recp_match = split(",", $all_recp_matches[0]);
        		if ($best_recp_match[1] eq $best_match[0]){
                    # se coge la menos restrictiva
        			$long_cut = $seqobj->length <= $seq_obj_look->length ? $seqobj->length : $seq_obj_look->length;
                } else{
                    next
                }
            } else{
                next
            }
        } else{
            # organismo no proporcionado, la fiabilidad del algoritmo se empobrece.
            $long_cut = $seqobj->length;
        }
		if ($best_match[2] >= 0.66*$long_cut  and $best_match[3] >= 66){
			# Alineamiento al menos del 66% de la mayor de las secuencias y similaridad de, al menos,
			# 66% -> tenemos un ortólogo (probablemente).
			$cont ++;
			my $self_id = $best_match[0];
			my $id_match = $best_match[1];
            &add_output($self_id, $id_match, "BDBH");
	    }
	}
    print "$cont matches where produced by BEST RECIPROCAL MATCHES\n";
    # Finalmente, se reinician los punteros de las streams
    foreach (@_) {
        seek($_->_fh, 0, 0);
    }
}

#4. Blast+
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
   system("$prog", "-db", "$path_now/$db_file", "-query", "$path_now/querytmp", "-max_target_seqs", "1", "-outfmt", "10 qseqid sseqid length ppos", "-out","$path_now/outmp");
   return "$path_now/outmp";
   # Re: eliminar archivos temporales en mainprint "$long_cut\n"; if ($i>20){last};
}

sub get_num_cores{
   # Función que busca el número de cores.
   # OUTPUT -> int, número de núcleos de procesador  o, si falla, 1.
   my $filename = '/proc/cpuinfo'; # archivo de información de linux
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

# Línea para hacer el script tanto módulo como ejecutable ("python-like")
&main() if not caller();

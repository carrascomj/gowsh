use warnings;
use strict;
use LWP::Simple;
use Getopt::Long;
use Data::Dumper;
use JSON;

our @path_tmp; # path a archivos temportales que serán borrados
our @path_files; # [0] gfile, [1] tfile y [2] modfile (de haberlo)
our $out_file; # almacena archivo de output
our %out_table; # almacena el output a lo largo del programa
our $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

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
    my $seqcount = `grep -c ">" $out`
    my $ret_point = 10000
    while (($seqcount / 10000) !~ /\./){ # while num secuencias es múltiplo de 10000
        my $new_url = $url . "&retstart=$retpoint";
        $out = get($new_url);
        open (OUT, ">>$name");
        print OUT $out;
        close OUT;
        $seqcount = `grep -c ">" $out`;
        $ret_point += 10000;
    }
    return $name;
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
    if (!$req_xml){
        return 0
    }
    # Parseamos el req_xml para buscar el(los) id(s) del genoma.
    my @UID;
    while ($req_xml =~ /<Id>(\d+?)<\/Id>/gs) {
        push @UID, $1
    }
    my $ids = join(',', @UID);
}

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

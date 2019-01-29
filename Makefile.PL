use 5.026001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'WebAPIsGOWSH',
    VERSION_FROM      => 'lib/WebAPIsGOWSH.pm', # finds $VERSION, requires EU::MM from perl >= 5.5
    PREREQ_PM         => { JSON => '2.97001', Bio::SeqIO => 1.007002}, # e.g., Module::Name => 1.1
    ABSTRACT_FROM     => 'lib/WebAPIsGOWSH.pm', # retrieve abstract from module
    AUTHOR            => 'Jorge <carrascomurielj@gmail.com>',
    #LICENSE           => 'MIT',
    #Value must be from legacy list of licenses here
    #http://search.cpan.org/perldoc?Module%3A%3ABuild%3A%3AAPI
);

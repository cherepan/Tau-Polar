#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$UserDir="";

if($UserID eq "vcherepa"){
    $UserIDCern="cherepan";
#    $UserDir="--vcherepa";
}

$PWD=getcwd;

printf("\n ---> Your user ID is:   $UserID \n");
if($ARGV[0] eq "--help" || $ARGV[0] eq ""){

    printf("\nThe installation follows the recomendation given in the description");
    printf("\nof LLRHiggsTauTau for the latest CMSSW release: https://github.com/LLRCMS/LLRHiggsTauTau");
    printf("\n\nThis code requires one input option. The syntax is: ./todo.pl [OPTION]");
    printf("\n\nRun todo script first with the '--settings <tauoladir>' optins to set up environment variable");
    printf("\n\nAfter this step is complet prcoceed further and ");
    printf("\npchoose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message\n");
    printf("\n./todo.pl --tauola  <tauoladir>           Install tauola \n");
    exit(0);  
}
my $dir = getcwd;
$time= strftime("%h_%d_%Y",localtime);



for($l=0;$l<$numArgs; $l++){
    
    if($ARGV[$l] eq "--settings"){
	$setdir=$ARGV[l+1];

	system(sprintf("rm Install_TauolaEnvironment_*"));	
	system(sprintf("echo \"export PYTHIA8DATA='$PWD/$setdir/tauola++/1.1.5/pythia8/176/xmldoc'\">> Install_TauolaEnvironment_$time"));
#	system(sprintf("echo \"export PYTHIA8DATA='/home-pbs/vcherepa/taua1/Installation/Tau-Polar/taudir/tauola++/1.1.5/pythia8/176/xmldoc'\" >> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PWD/TauSpiner/tauola++/1.1.5/lib\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PWD/$setdir/tauola++/1.1.5/pythia8/176/lib/\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PWD/$setdir/tauola++/1.1.5/HepMC-2.06.05/workdir/lib \">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PWD/$setdir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PWD/$setdir/tauola++/1.1.5/examples/UserCode\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"cernlib-use --version 5.34.18 root \n\">> Install_TauolaEnvironment_$time"));

  

	printf("\n\nInstructions:");
	printf("\nTo complete this step do source Install_TauolaEnvironment_$time \n\n");
    }
    if($ARGV[$l] eq "--tauola"){
	$tauoladir=$ARGV[l+1];
#	$l++;
	$currentdir=getcwd;


	system(sprintf("cernlib-use --version 5.34.18 root \n"));
                                

	system(sprintf("rm -rf $tauoladir \n"));
	printf("\nInstalling Tauola++  to  $tauoladir \n");
	system(sprintf("mkdir  $tauoladir \n"));
	system(sprintf("cd $tauoladir; wget http://service-spi.web.cern.ch/service-spi/external/MCGenerators/distribution/tauola++/tauola++-1.1.5-src.tgz; tar -xzvf tauola++-1.1.5-src.tgz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://service-spi.web.cern.ch/service-spi/external/MCGenerators/distribution/pythia8/pythia8-176-src.tgz ; tar -xzvf pythia8-176-src.tgz ;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://www.hepforge.org/archive/lhapdf/lhapdf-5.9.1.tar.gz; tar -xzvf lhapdf-5.9.1.tar.gz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://mc-tester.web.cern.ch/MC-TESTER/MC-TESTER-1.25.0.tar.gz; tar -xzvf MC-TESTER-1.25.0.tar.gz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://lcgapp.cern.ch/project/simu/HepMC/download/HepMC-2.06.05.tar.gz;  tar -xzvf  HepMC-2.06.05.tar.gz;"));
	printf("\n Downloading Complete ... \n");
	printf("\n Start Installation ... \n\n\n");
	printf("\n ___________________Installing HEPMC ... _____________________\n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; .././configure -prefix=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir  -with-momentum=GEV -with-length=CM; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; make; make check; make install;"));
	printf("\n___________________Installing  LHAPDF ... _____________________\n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/; ./configure -prefix=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir  --libdir=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/; make;  make install;"));
	system(sprintf("mkdir  $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/share/lhapdf/PDFsets;"));
	system(sprintf("cd  $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/share/lhapdf/PDFsets; ../../../bin/lhapdf-getdata --repo=http://www.hepforge.org/archive/lhapdf/pdfsets/5.9.1 MSTW2008nnlo90cl.LHgrid; "));
	printf("\n ___________________Installing pythia8 ... _____________________ \n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/;  ./configure --enable-shared --lcgplatform=slc6_amd64_gcc530-opt --with-hepmc=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir --with-hepmcversion=2.06.05;"));	
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/;  make; "));
	printf("\n ___________________Compiling MC-TESTER ... _____________________\n\n\n");
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/;  ./configure --with-HepMC=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/;  make; "));
	printf("\n ___________________Compiling tauola ... _____________________\n\n\n");

	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export PYTHIA8DATA=$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/xmldoc;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauSpiner/tauola++/1.1.5/lib;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/lib/;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/;export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir/lib; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib"));
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/;   ./configure --prefix=$PWD/$tauoladir/tauola++/1.1.5/workdir  --with-hepmc=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir  --with-pythia8=$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/  --with-lhapdf=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/ --with-mc-tester=$PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/   --with-tau-spinner; "));
	system(sprintf("cp Makefile  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
	system(sprintf("cp mypythia_example.cxx  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
	system(sprintf("echo \"\#dummy\">> $PWD/$tauoladir/tauola++/1.1.5/examples/.deps/mypythia_example.Po"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples;  mkdir UserCodes;"));
	system(sprintf("cp UserCodes/Makefile.shared $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/MultiplyNumbers.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/MultiplyNumbers.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/AddNumbers.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/AddNumbers.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/main.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes/; make all -f Makefile.shared;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make all;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make install;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples; make"));
    }

}


#!/usr/bin/perl
#

# Do svn operations to produce blame files

#  while (glob("*.f90 *.c *.h"))
#    {
#      $infilename = $_;
#      print "infilename = $infilename \n";
#      @fileparts = split /.f90/, $infilename;
#      $outfilename = @fileparts[0].".svnblame";
#      $command = "rm -rf $outfilename";
#      system( $command );
#      $command = "svn blame $infilename > $outfilename";
#      print "command = $command \n";
#      system( $command );
#    }
#

# Parse each .svnblame file

%developerlines = ();
while (glob("*.svnblame"))
  {
   $infilename = $_;
   print "infilename = $infilename \n";
   open(infile,"<$infilename");
   while ( $infileline = <infile> )
     {
       chomp($infileline);
       @infilelineparts = split /\s+/, $infileline;
       $developer = $infilelineparts[2];
       $developerlines{$developer}++;
     }
   close(infile);
  }

# Set up developer information

%developerinfo = (
  "abowery" => {"name" => "Andy Bowery", "institution" => "Oxford/KCL"},
  "adv81" => {"name" => "Adelaide de Vecchi", "institution" => "Oxford/KCL"},
  "ahungnz" => {"name" => "Alice Hung", "institution" => "ABI"},
  "andync42" => {"name" => "Andrew Cookson", "institution" => "Oxford/KCL"},
  "areeve" => {"name" => "Adam Reeve", "institution" => "ABI"},
  "catalept" => {"name" => "Caton Little", "institution" => "ABI"},
  "chloeiw" => {"name" => "Chloe Irwin Whitney", "institution" => "ABI"},
  "chrispbradley" => {name => "Chris Bradley", "institution" => "ABI"},
  "chrm76" => {"name" => "Christian Michler", "institution" => "Oxford/KCL"},
  "cistib" => {"name" => "Ali Pashaei", "institution" => "UPF"},
  "dbrooks12" => {"name" => "Dave Brooks", "institution" => "ABI" },
  "dougalcowan" => {"name" => "Dougal Cowan", "institution" => "ABI" },
  "eltomo" => {"name" => "Thomas Heidlauf", "institution" => "Stuttgart" },
  "ir10"=> {"name" => "Ishani Roy", "institution" => "Oxford/KCL" },
  "isotopeb12" => {"name" => "Bojan Blazevic", "institution" => "Oxford/KCL" },
  "jackisawesome" => {"name" => "Jack Lee", "institution" => "Oxford/KCL" },
  "jagirhussan" => {"name" => "Jagir Hussan", "institution" => "ABI" },
  "jbudelmann" => {"name" => "Jin Budelmann", "institution" => "ABI" },
  "kmith" => {"name" => "Kuman Mithraratne", "institution" => "ABI" },
  "krittian" => {"name" => "Sebastian Krittian", "institution" => "Oxford/KCL" },
  "marktrew" => {"name" => "Mark Trew", "institution" => "ABI" },
  "martinsteg" => {"name" => "Martin Steghofer", "institution" => "UPF" },
  "mmcc094" => {"name" => "Matthew McCormick", "institution" => "Oxford/KCL" },
  "nickerso" => {"name" => "David Nickerson", "institution" => "ABI" },
  "nordslet23" => {"name" => "David Nordsletten", "institution" => "Oxford/KCL" },
  "olewsaa" => {"name" => "Ole Saastad", "institution" => "Oslo" },
  "onordb" => {"name" => "Oyvind Nordbo", "institution" => "UMB" },
  "oroh001" => {"name" => "Oliver Rohrle", "institution" => "Stuttgart" },
  "peterhunter" => {"name" => "Peter Hunter", "institution" => "ABI" },
  "prasad012" => {"name" => "Prasad Gamage", "institution" => "ABI" },
  "robertjiahe" => {"name" => "Robert Jiahe", "institution" => "Oxford/KCL" },
  "sanderland" => {"name" => "Sander Land", "institution" => "Oxford/KCL" },
  "snie007" => {"name" => "Steven Niederer", "institution" => "Oxford/KCL" },
  "soroushsafaei" => {"name" => "Soroush Safaei", "institution" => "ABI" },
  "tingy" => {"name" => "Ting Yu", "institution" => "ABI" },
  "ttme" => {"name" => "Tim Wu", "institution" => "ABI" },
  "vizaxis" => {"name" => "Randall Brittan", "institution" => "ABI" },
  "wolfye" => {"name" => "Heye Zhang", "institution" => "ABI" },
  "xlu012" => {"name" => "Xiao Bo Lu", "institution" => "ABI" }
 );

# Find the total number of lines and the number of lines for each institution

$totalnumberoflines = 0;
%institutionlines = ();
foreach $developer (keys %developerlines)
  {
    $totalnumberoflines += $developerlines{$developer};
    $institution = $developerinfo{$developer}{"institution"};
    $institutionlines{$institution} += $developerlines{$developer}
  }

# Output the summary information.

open(outfile,">developerslines.txt");

printf outfile "The total number of lines is = %7d\n",$totalnumberoflines;
printf outfile "\n";

printf outfile "%15s     %20s     %15s  %15s\n","Sourceforge ID","Name","# of Lines","Percentage";
printf outfile "\n";
foreach $developer (keys %developerlines)
  {
    $percentage = 100.0*$developerlines{$developer}/$totalnumberoflines;
    $name = $developerinfo{$developer}{"name"};
    printf outfile "%15s     %20s     %15d        %7.2f %\n",$developer,$name,$developerlines{$developer},$percentage;
  }
printf outfile "\n";
printf outfile "%15s                              %15s  %15s\n","Institution","# of Lines","Percentage";
printf outfile "\n";
foreach $institution (keys %institutionlines)
  {
    $percentage = 100.0*$institutionlines{$institution}/$totalnumberoflines;
    printf outfile "%15s                              %15d        %7.2f %\n",$institution,$institutionlines{$institution},$percentage;
  }

close(outfile);

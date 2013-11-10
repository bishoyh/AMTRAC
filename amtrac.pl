#!/usr/local/bin/perl
# AMTRAC -- AMazing TRimming And Cleaning ---
# Use to trim barcodes and other motifs from DNA sequences
# Usage:
# amtrac.pl -b barcodes.fasta fastafile
# amtrach.pl -a GCCCC -b barcodes.fasta fastafile # -a can specify a changing anchor
#
# Copyright (C) 2013 Bishoy Hanna, Elizabeth Green and Monica Medina
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.



use Term::ANSIColor;
use Getopt::Long;
print STDERR colored( "------AMTRAC-------\n AMazing TRimming And Cleaning", 'bright_red' ), "\n";
$heredoc = <<END;

       (@@@)     (@@@@@)
                               (@@)     (@@@@@@@)        (@@@@@@@)
                         (@@@@@@@)   (@@@@@)       (@@@@@@@@@@@)
                    (@@@)     (@@@@@@@)   (@@@@@@)             (@@@)
               (@@@@@@)    (@@@@@@)                (@)
           (@@@)  (@@@@)           (@@)
        (@@)              (@@@)
       .-.               
       ] [    .-.      _    .-----.
     ."   """"   """""" """"| .--`
    (:--:--:--:--:--:--:--:-| [___    .------------------------.
     |C&O  :  :  :  :  :  : [_9_] |'='|.----------------------.|
    /|.___________________________|___|'--.___.--.___.--.___.-'| 
   / ||_.--.______.--.______.--._ |---\'--\-.-/==\-.-/==\-.-/-'/--
  /__;^=(==)======(==)======(==)=^~^^^ ^^^^(-)^^^^(-)^^^^(-)^^^ jgs
~~~^~~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~~^~~
END
 print STDERR colored($heredoc, bright_blue);
 
$ANCHOR = 'CACGGAG';			# Default Anchor something down stream if needed, the script can still just use everything after the barcodes
$::BARCODES_FILE = '';		# -f, name of file containing a list of barcodes in fasta format
my %opts = (config => '');

GetOptions(\%opts, qw(
   barcodesfile|b=s
   anchor|a
));

$::BARCODES_FILE=$opts{barcodesfile};
use Data::Dumper; # for debugging arrays
  open(BARCODES, $::BARCODES_FILE) or die("Cannot open BARCODES_FILE '$::BARCODES_FILE', $!\n");
  print STDERR "BARCODESFILE='$::BARCODES_FILE'\n" ;
  while (defined ($line = <BARCODES>))
    {
      

    chomp $line;

        if ($line =~ /^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $line;
            $seqData[$i] = "";
        } else {
            s/^\s+//; s/\s+$_//;
	    s/\s+//g;                  # get rid of any spaces
            next if (/^$line/);            # skip empty line
            $seqData[$i] = $seqData[$i] . uc($line);
        }

push (@BarCodesList, $line);
if ($line =~ /^>/){
  print STDERR colored("Barcode ID = $line\t", bright_cyan);
}else{
   print STDERR colored("Sequence ='$line'\n" , bright_magenta);
    }
 }

    
  close BARCODES;

### start reading the fasta file -- main program loop
$fasta_input_file = shift @ARGV;
$fasta_input_file = '-' if (!$fasta_input_file);
open(FASTAIN, $fasta_input_file) || die("Can't open fasta_input_file: '$fasta_input_file'\n");
  chomp $line;  chomp $line;
while ($line = <FASTAIN>)
  {
  chomp $line;
  $line_num++;
  if ($line =~ /^>/)
    {
    if ($found_data >= 0)	# check if we are running for the first time - parsing the fasta record
      {
     trim_contig($contig, $sequence) if ($contig);
      }
    $found_data = 0;		# reset the counter
   
    $header = $line;
    if ($header =~ /^>((\S+).*)/)
      {
      
        $contig = $1;
       
      }
    else
      {
      print STDERR "Error: Invalid fasta_input_file format: '$fasta_input_file'\n  Fasta input line number=$line_num\n";
      exit 2;
      }
    $sequence = '';
    }
  else 
    {
    if ($found_data < 0)
      {
      print STDERR "Error: Invalid fasta_input_file format: '$fasta_input_file'\n";
      exit 2;
      }
    $line =~ s/\s//g;
    $sequence .= $line;
    $found_data = 1;
    }
  } # end while
if ($contig)
  {
  trim_contig($contig, $sequence);
  }
else
  {
  print STDERR "Error: Empty fasta_input_file: '$fasta_input_file'\n";
  exit 2;
  }
print STDERR "\nEnd of searches\n";

close(FASTAIN);


### Do the actual matching and trimming

sub trim_contig
  {
  my($contig, $sequence) = @_;
  my($pattern, $len, $fhits, $rhits, $barcode_start, $barcode_end, $barcode_len, $match);
  $len = length($sequence);
  $fhits = $rhits = 0;

  print STDERR colored( "Searching $len bases of\t$contig\n" , bright_yellow);

  my $num_patterns = scalar @BarCodesList;
  foreach $pattern (@BarCodesList)
    {
if ($pattern =~ /^>/){
	$barcodeID=$pattern;
next;
	}
    while ($sequence =~ /($pattern)/ig)
      {

@Anchorindex = match_all_positions($ANCHOR, $sequence);
@Barcodeindex = match_all_positions($pattern, $sequence);
print STDERR colored("##### Barcode Ends here $Barcodeindex[0][1]\t Anchor starts here $Anchorindex[0][0]\n",bright_green);
      $match = $1;
	$seqlen= length $sequence;
	$anchoneg = $seqlen - $Anchorindex[0][0];
      $barcode_len = length $match;
      $barcode_end = pos $sequence;
      $barcode_start = $barcode_end - $barcode_len + 1;
      pos $sequence -= ($hit_len - 1)
        if ($::ALLOW_OVERLAPS && $barcode_len > 1);
	if ($Anchorindex[0][0] > 1){
		
          print "$barcodeID\_$contig\n" , substr($sequence,$barcode_end,-$anchoneg),"\n";
          $fastarec = "$barcodeID\_$contig\n" . substr($sequence,$barcode_end,-$anchoneg) ."\n";
		
                append_to_fasta($barcodeID,$fastarec);
	}else{
	print STDERR colored("##########Anchor not found#######\n", bright_magenta);
      print "$barcodeID\_$contig\_NoAnchor\n" , substr($sequence,$barcode_end),"\n";
	$fastarec = "$barcodeID"."_"."$contig"."_NoAnchor\n".substr($sequence,$barcode_end)."\n";
append_to_fasta($barcodeID,$fastarec);

}
$fhits++;
      }

    }
  } # end trim_contig


#############M#####################################
#	Match All Positions SubRoutine            #
###################################################

sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, [ $-[0], $+[0] ];
    }
    return @ret
}

########Subroutine to append to fastafiles########
sub append_to_fasta {

my ($filename,$fastarecord) = @_;
$filename=~ tr/>//d;
open(my $fh, '>>', $filename) or die "Could not open file '$filename' $!";
print $fh $fastarecord;
close $fh;

}



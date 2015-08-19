#!/bin/perl
#perl parseFBvcf.pl /home/tng/proj/Battle_GC/case4/run*/analysis/geno/freebayes/unfiltered.vcf |less

use warnings;
use strict;

  local $SIG{__WARN__} = sub {
    my $message = shift;
    logger('warning', $message);
    logger('in', $_);
  };
 
  my $counter;
  count();
  print "$counter\n";
  sub count {
    $counter = $counter + 42;
  }
 
 
  sub logger {
    my ($level, $msg) = @_;
    if (open my $out, '>>', 'log.txt') {
        chomp $msg;
        print $out "$level - $msg\n";
    }
  }

while(<>) {
	next if /^#/;
	my @l = split "\t"; 
	my @ancAllele = split "", $l[3];
	my @derivAllele = split ",", $l[4];
        my $pos = $l[1];
	my $contig = $1 if $l[0]=~/contig_(\d*)/;

        my $allAllele;
	push @{$allAllele}, \@ancAllele;
        my $len = $#ancAllele ;
	my $sameLen = 0;
	for (my $i = 0; $i <= $#derivAllele ; $i++) {
		my @temAllele = split "", $derivAllele[$i];
		push @{$allAllele}, \@temAllele;
		$samLen = 1 if $len != $#temAllele;
                $len = $#temAllele if $len < $#temAllele; 
	}
	
	if ($len > 0) {
		if ($samLen == 0) {
			# handling case just for pure haplotype subsitution (w/o any del or sub) 
			for (my $i = 0; $i <= $len; $i++) {
				my $isSNP = 0;
                                my $refNucl = $allAllele->[0][$i];
				my $devNucl="err";
				for (my $j = 1; $j <= $#{$allAllele}; $j++) {
					($devNucl, $isSNP) = ($allAllele->[$j][$i], 1) if $refNucl ne $allAllele->[$j][$i];	
				}
				if ($isSNP == 1) {
					print join "\t", $contig, $pos+$i, $refNucl, $devNucl; 
					map{print "\t", $allAllele->[$1][$i], $allAllele->[$2][$i]  if $l[$_] =~/^([01])\/([01])/}(9..$#l);
                			print "\n";
				}
			}
		}
		else {
			if ($l[7]=~/del/) {
				#deletion case: minus 1 bp
				my $lenAlt = $#derivAllele + 1;
                                my @nucl = ($allAllele->[0][1], (".")x$lenAlt) ;
				for(my $j = 1; $j <= $#{$allAllele}; $j++) {
					$nucl[$j] = $allAllele->[0][1] if $#{$allAllele->[$j]} == $#{$allAllele->[0]};
				}
				print join "\t", $contig, $pos+1, $nucl[0], "."; 
				map{print "\t", $nucl[$1], $nucl[$2]  if $l[$_] =~/^([01])\/([01])/}(9..$#l);
                		print "\n";

				for (my $i = 2; $i <= $len; $i++) {
					my $isSNP = 0;
                                	my $refNucl = $allAllele->[0][$i];
					my $devNucl;
					my @nucl;
					push @nucl, $refNucl;
					for (my $j = 1; $j <= $#{$allAllele}; $j++) {
						my $adj = $i - ($#{$allAllele[0]} - $#{$allAllele[$j]});
						push @nucl, $allAllele->[$j][$adj]
						($devNucl, $isSNP) = ($allAllele->[$j][$adj], 1) if $refNucl ne $allAllele->[$j][$adj];	
					}
					if ($isSNP == 1) {
						print join "\t", $contig, $pos+$i, $refNucl, $devNucl; 
						map{print "\t", $nucl[$1], $nucl[$2]  if $l[$_] =~/^([01])\/([01])/}(9..$#l);
                				print "\n";
					}
				}
			}
			elsif($l[7]=~/ins/){
				print $_, "\n";
			}
			else{
				print $_,"\n";
			}
		}

	}
	else{
		print join "\t", $contig, $pos, @{$allAllele->[0]}, @{$allAllele->[1]}; 
		map{print "\t", @{$allAllele->[$1]}, @{$allAllele->[$2]}  if $l[$_] =~/^([01])\/([01])/}(9..$#l);
                print "\n";
	}
}

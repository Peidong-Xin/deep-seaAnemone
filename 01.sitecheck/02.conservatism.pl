use strict;
use warnings;

## Conservativeness of sites assessed by Shannon entropy H(x)=-sum[P(x)*log2P(x)]
## Placeholder of gap (-) is also considered to be an AA state. But if the number of gaps exceeds 60% of the total, this site will no longer be calculated and assigned the value x.
## Subjective interpretation: A value of 0 means extremely conservative, less than 0.5 can be considered conservative, less than or equal to 1 is considered tolerable, greater than 1 is not conservative, and x is not considered.

my $dir="01.alignment";
my @file=<$dir/*.fa>;
my $sum_file=@file;
my $do_file=0;
open (O, "> $0.sev");

foreach my $file(@file){
	$do_file++;
    my %AA;
    open (I, "$file");
    $/ = ">";
    my $len;
    my @sum_sp;
    
    while(<I>){
	chomp;
	next if (/^$/);
	my @a=split(/\n/, $_);
	my $titou=shift @a;
	$titou=~s/^\s+|\s+$//g;
	my ($id)=split(/\s+/, $titou);
	die "Unable to parse a sequence identifier in $file\n" unless defined $id and length $id;
	push @sum_sp, $id;
	my $seq=join "",@a;
	my @base=split(//,$seq);
	$len=@base;
	for(my $i=0;$i<$len;$i++){
	    $AA{$i}{$base[$i]}++;
	}
    }
    close I;
    $/ = "\n";
    
    my $sum=@sum_sp;
    my (@site, @shannon);
    my $line=$sum * 0.6;
    my $gap="-";
    
    for(my $i=0; $i<$len; $i++){
	my $order=$i+1;
	push @site, $order;
	if($AA{$i}{$gap}){
	    if($AA{$i}{$gap} > $line){
		push @shannon, "x";
	    }
	    else{
		my $Hx=0;
		foreach my $aa(keys %{$AA{$i}}){
		    my $rate=$AA{$i}{$aa} / $sum;
		    my $Px= log($rate) / log(2);
		    $Hx += -($rate * $Px);
		}
		$Hx=sprintf"%.3f",$Hx;
		push @shannon, $Hx;
	    }
	}
	else{
	    my $Hx=0;
	    foreach my $aa(keys %{$AA{$i}}){
		my $rate=$AA{$i}{$aa} / $sum;
		my $Px= log($rate) / log(2);
		$Hx += -($rate * $Px);
	    }
	    $Hx=sprintf"%.3f",$Hx;
	    push @shannon, $Hx;
	}
    }
    
    print O "## $file (length $len)\n";
    print O join "\t",@site,"\n";
    print O join "\t",@shannon,"\n";
    print O "\n";
	## In-situ display progress
	local $| = 1;
	my $bili=($do_file/$sum_file)*100;
	$bili=sprintf"%.2f",$bili;
	print "\r $bili"."%"."completed";
	local $| = 0;
}
close O;

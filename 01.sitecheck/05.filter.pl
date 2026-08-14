use strict;
use warnings;


my @files=<"01.alignment/*.fa">;
my $info="./Nevc.lst"; # Target-gene annotation information
my @infs=<"04.group_and_Theta-w.pl.20.strict">;
my $thre=0.6; ## Theta-w greater than $thre was picked

my (%hash, %have);

open I, "$info";
while(<I>){
	chomp;
	next if(/^#/);
	my @a=split(/\s+/);
	$a[2]=~/^([^-]+)/;
	$hash{$a[1]}="$a[1]\t$1";
}
close I;

my $time1=localtime();
print "$time1 gene info done!\n";

foreach my $file(@files){
	open I, "$file";
	while(<I>){
		chomp;
		next unless (/^>/);
		$_=~/>(\S+)/;
		if($hash{$1}){
			$have{$file}=$hash{$1};
		}
	}
	close I;
}

my $time2=localtime();
print "$time2 gene assign to fas file!\n";

foreach my $inf(@infs){
	$inf=~/(\w+)$/;
	my $type=$1;
	my %harv;
	open I, "$inf";
	open O, "> $0.$thre.$type";
	$/="##";
	while(<I>){
		chomp;
		next if(/^$/);
		my @a=split(/\n/);
		my $head=shift @a;
		my $text=join "\n", @a;
		$head=~/(\w+)\s+([\.0-9]+)/;
		my $group=$1;
		my $val=$2;
		if($val > $thre){
			next unless($text=~/^(.+?\.fa)\s+\(Site\s+\d+\):/m);
			my $target_file=$1;
			my $nt="## $group\t$have{$target_file}\t$val\n$text";
			$harv{$val}{$group}=$nt;
		}
	}
	close I;
	$/="\n";
	foreach my $k(sort {$b <=> $a} keys %harv){
		foreach my $g(sort {$a cmp $b} keys %{$harv{$k}}){
			print O "$harv{$k}{$g}\n\n";
		}
	}
	close O;
	print "$inf was done!\n";
}

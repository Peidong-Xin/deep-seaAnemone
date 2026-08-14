use strict;
use warnings;

my $shre=shift;
$shre=10 unless defined $shre;
die "$0 <window size; use 10 for the STAR Methods analysis>\n" unless $shre=~/^\d+$/ and $shre>0;
my @files=<"03.detect_site_variation.pl.strict">;
my $shannon="02.conservatism.pl.sev";
my $thd=5; # Loci less than $thd were evaluated as a whole

## Theta-w = (AA value)/(counted AA number), AA value, viz. each site was assigned a value based on the shannon index (SI): SI(0)=1, SI(<=0.5)=0.5, SI(<=1)=0.3.

my $statLen=$shre * 2;
my (%hash, %have);

foreach my $file(@files){
	$file=~/(\w+)$/;
	my $type=$1;
	my $group=0;
	open I, "$file";
	while(<I>){
		chomp;
		next unless(/^(.+?\.fa)\s+\(Site\s+(\d+)\):/);
		my ($gen, $site)=($1, $2);
		$have{$gen}=1;
		if($hash{$type}{$group}){
			if($hash{$type}{$group}{$gen}){
				my @weizhi = sort {$b <=> $a} keys %{$hash{$type}{$group}{$gen}};
				my $far=abs($site-$weizhi[0]-1);
				#print "$_ : $gen | $site | $far\n";
				if($far < $thd){
					$hash{$type}{$group}{$gen}{$site}=$_;
				}
				else{
					$group++;
					$hash{$type}{$group}{$gen}{$site}=$_;
				}
			}
			else{
				$group++;
				$hash{$type}{$group}{$gen}{$site}=$_;
			}
		}
		else{
			$group++;
			$hash{$type}{$group}{$gen}{$site}=$_;
		}
	}
	close I;
}

##Debug
=pod
foreach my $type(keys %hash){
	open O, "> $0.test.$type";
	foreach my $group(sort{$a <=> $b} keys %{$hash{$type}}){
		print O "##group$group\n";
		foreach my $gen(keys %{$hash{$type}{$group}}){
			foreach my $site(sort{$a <=> $b} keys %{$hash{$type}{$group}{$gen}}){
				print O "$hash{$type}{$group}{$gen}{$site}\n";
			}
		}
	}
	close O;
}
=cut

print "group and merge were done | and now begin $shannon\n";

my %SEI;
open I, "$shannon";
while(<I>){
	chomp;
	next unless(/^##\s+(.+?\.fa)\s+\(length\s+\d+\)/);
	my $gen=$1;
	next unless($have{$gen});
	my $loc=<I>;
	my $val=<I>;
	chomp ($loc, $val);
	my @sites=split(/\s+/, $loc);
	my @values=split(/\s+/, $val);
	my $len=@sites;
	for(my $i=0;$i<$len;$i++){
		$SEI{$gen}{$sites[$i]}=$values[$i];
	}
}
close I;

print "The processing of the $shannon has been completed. The theta count now begins\n";

foreach my $type(keys %hash){
	open O, "> $0.$statLen.$type";
	foreach my $group(sort {$a <=> $b} keys %{$hash{$type}}){
		my @genes=keys %{$hash{$type}{$group}};
		print "ERROR: more than one genes were found in a group, DELETE the results and CHECK perl again!!!\n" if (@genes > 1);
		my $gen=shift @genes;
		die "No conservation data were found for $gen\n" unless exists $SEI{$gen};
		my @locs= sort {$a <=> $b} keys %{$hash{$type}{$group}{$gen}};
		my @sts=@locs;
		if(@locs > 1){
			my $t_min=shift @locs;
			my $t_max=pop @locs;
			my $up=$t_min-$shre;
			my $down=$t_max+$shre;
			my $sum=0;
			my $stat=0;
			if($up < 1){
				$up=1;
			}
			my @bases=keys %{$SEI{$gen}};
			if($down > @bases){
				$down=@bases;
			}
			for(my $i=$up;$i<$t_min;$i++){
				next unless exists $SEI{$gen}{$i};
				next if($SEI{$gen}{$i} eq "x");
				$sum++;
				$stat+=1 if($SEI{$gen}{$i} < 0.01);
				$stat+=0.5 if($SEI{$gen}{$i} <= 0.5 and $SEI{$gen}{$i} >= 0.01);
				$stat+=0.3 if($SEI{$gen}{$i} <= 1 and $SEI{$gen}{$i} > 0.5);
			}
			for(my $i=$t_max+1;$i<=$down;$i++){
				next unless exists $SEI{$gen}{$i};
				next if($SEI{$gen}{$i} eq "x");
				$sum++;
				$stat+=1 if($SEI{$gen}{$i} < 0.01);
				$stat+=0.5 if($SEI{$gen}{$i} <= 0.5 and $SEI{$gen}{$i} >= 0.01);
				$stat+=0.3 if($SEI{$gen}{$i} <= 1 and $SEI{$gen}{$i} > 0.5);
			}
			die "No eligible flanking sites were found for $gen group$group\n" if $sum==0;
			my $theta=$stat/$sum;
			$theta=sprintf "%.3f",$theta;
			print O "## group$group\t$theta\n";
			foreach my $lo(@sts){
				print O "$hash{$type}{$group}{$gen}{$lo}\n";
			}
			my @weizhi=($up..$down);
			my @rate;
			foreach my $k(@weizhi){
				push @rate, $SEI{$gen}{$k};
			}
			print O join "\t",@weizhi,"\n";
			print O join "\t",@rate,"\n\n";
		}
		else{
			my $targ=shift @locs;
			my $up=$targ-$shre;
			my $down=$targ+$shre;
			my $sum=0;
			my $stat=0;
			if($up < 1){
				$up=1;
			}
			my @bases=keys %{$SEI{$gen}};
			if($down > @bases){
				$down=@bases;
			}
			for(my $i=$up;$i<$targ;$i++){
				next unless exists $SEI{$gen}{$i};
				next if($SEI{$gen}{$i} eq "x");
				$sum++;
				$stat+=1 if($SEI{$gen}{$i} < 0.01);
				$stat+=0.5 if($SEI{$gen}{$i} <= 0.5 and $SEI{$gen}{$i} >= 0.01);
				$stat+=0.3 if($SEI{$gen}{$i} <= 1 and $SEI{$gen}{$i} > 0.5);
			}
			for(my $i=$targ+1;$i<=$down;$i++){
				next unless exists $SEI{$gen}{$i};
				next if($SEI{$gen}{$i} eq "x");
				$sum++;
				$stat+=1 if($SEI{$gen}{$i} < 0.01);
				$stat+=0.5 if($SEI{$gen}{$i} <= 0.5 and $SEI{$gen}{$i} >= 0.01);
				$stat+=0.3 if($SEI{$gen}{$i} <= 1 and $SEI{$gen}{$i} > 0.5);
			}
			die "No eligible flanking sites were found for $gen group$group\n" if $sum==0;
			my $theta=$stat/$sum;
			$theta=sprintf "%.3f",$theta;
			print O "## group$group\t$theta\n$hash{$type}{$group}{$gen}{$targ}\n";
			my @weizhi=($up..$down);
			my @rate;
			foreach my $k(@weizhi){
				push @rate, $SEI{$gen}{$k};
			}
			print O join "\t",@weizhi,"\n";
			print O join "\t",@rate,"\n\n";
		}
	}
	close O;
	print "$0.$statLen.$type was done\n";
}

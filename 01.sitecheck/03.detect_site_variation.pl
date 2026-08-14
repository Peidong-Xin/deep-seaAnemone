use strict;
use warnings;

## modified 2019-6-11 and 2022-4-20
## STRICT: There are two kinds of amino acid state in total, target species occupies only one state while other species occupies another one
## CONSTRAINT: There are more than two kinds of amino acid state in total, target species occupies only one state while other species occupies others but does not share with target species
## LOOSE: There are more than two kinds of amino acid state in total, target species occupies more than one state while other species occupies another one that is not shared with target species

my $dir="01.alignment";
my @file=<$dir/*.fa>;
my %ingroup=('Target_deep_sea_sp1'=>'a', 'Target_deep_sea_sp2'=>'a', 'Target_deep_sea_sp3'=>'a', 'Target_deep_sea_sp4'=>'a');

sub is_ingroup {
    my ($id)=@_;
    foreach my $target (keys %ingroup){
        return 1 if $id eq $target or $id=~/^\Q$target\E_/;
    }
    return 0;
}

open (O, "> $0.strict");
open (L, "> $0.loose");
open (C, "> $0.constraint");
foreach my $file(@file){
my %hash;
open (I, "$file");
$/ = ">";
my $len;

while(<I>){
    chomp;
    next if (/^$/);
    my @a=split(/\n/, $_);
    my $titou=shift @a;
	$titou=~s/^\s+|\s+$//g;
	my ($id)=split(/\s+/, $titou);
	die "Unable to parse a sequence identifier in $file\n" unless defined $id and length $id;
    my $seq=join "",@a;
    my @base=split(//,$seq);
    $len=@base;
    #@{$hash{$id}}=@base;
    for(my $i=0;$i<$len;$i++){
        $hash{$i}{$id}=$base[$i];
    }
}
close I;
$/ = "\n";

for(my $i=0; $i<$len; $i++){
    my (%site, %test, %cov);
    # print "$i\n";
    foreach my $id(keys %{$hash{$i}}){
        my $base=$hash{$i}{$id};
        $site{$base}{$id}="yes";
        if (is_ingroup($id)){
            $cov{$base}{$id}="yes";
        }
        else{
            $test{$base}{$id}="yes";
        }
    }
    my @stat=keys %site;
    # print "$i\n";
    # print join "\t",@stat,"\n";

    if (@stat == 2){
        if(keys %test == 1){
	if(keys %cov == 1){
            my @aa=keys %test;
            my $aa=shift @aa;
			#next if ($aa eq "-");
            my @bb=keys %cov;
            my $bb=shift @bb;
			#next if ($bb eq "-");
            my $order=$i+1;
            my (@title, @nu);
            foreach my $k(sort {$a cmp $b} keys %{$hash{$i}}){
	push (@title, $k);
	push (@nu, $hash{$i}{$k});
            }
            print O "$file (Site $order): $aa ==> $bb\n";
            print O join "\t",@title,"\n";
            print O join "\t",@nu,"\n\n";
	}        
	else{
	next;
	}
	}
        else{
            next;
        }
    }
	elsif(@stat > 2){
		if(keys %test == 1){
			my @aa=keys %test;
			my $aa=shift @aa;
			next if ($aa eq "-");
			my %lnshi;
			foreach my $bb(keys %cov){
				$lnshi{$bb}=1;
			}
			next if($lnshi{$aa});
			my $order=$i+1;
			my (@title, @nu);
			foreach my $k(sort {$a cmp $b} keys %{$hash{$i}}){
				push (@title, $k);
				push (@nu, $hash{$i}{$k});
			}
			my $shul=@stat - 1;
			print L "$file (Site $order): $aa ==> $shul\n";
			print L join "\t",@title,"\n";
			print L join "\t",@nu,"\n\n";
		}
		else{
			next unless (keys %cov == 1);
			my @bb=keys %cov;
			my $bb=shift @bb;
			next if ($bb eq "-");
			my %lnshi;
			foreach my $aa(keys %test){
				$lnshi{$aa}=1;
			}
			next if($lnshi{$bb});
			my $order=$i+1;
			my (@title, @nu);
			foreach my $k(sort {$a cmp $b} keys %{$hash{$i}}){
				push (@title, $k);
				push (@nu, $hash{$i}{$k});
			}
			my $shul=@stat - 1;
			print C "$file (Site $order): $shul ==> $bb\n";
			print C join "\t",@title,"\n";
			print C join "\t",@nu,"\n\n";
		}
	}
    else{
        next;
    }
}
}
close O;
close L;
close C;

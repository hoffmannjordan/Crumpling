#!/usr/bin/perl
use Getopt::Std;

getopts("hg");

if(defined($opt_h)) {
	print "Usage: mesh2cyl <switches> <mesh2_file> <cylinder_file>\n\n";
	print "Switches:\n";
	print "-g                (Output in Gnuplot instead of POV format)\n";
	print "-h                (Print this information)\n";
	print "-u                (Add union keywords to POV output\n";
	exit 0;
}

die "Two arguments required" unless @ARGV==2;

$l=0;$vs[0]=0;
open A,"$ARGV[0]";
while(<A>) {
	next unless /mesh2/;
	$_=<A>;die "Incorrect header\n" unless /vertex_vectors/;
	$_=<A>;die "Incorrect header\n" unless /^\t\t(\d*),$/;

	# Read in the position vectors
	$n=$1;
	foreach $i ((1+$vs[$l])..($n+$vs[$l])) {
		$_=<A>;
		/^\t\t<(.*),(.*),(.*)>,?$/ or die "Unreadable position vector\n";
		$x[$i]=$1;
		$y[$i]=$2;
		$z[$i]=$3;
		$v[$i][0]=0;
	}

	# Skip the normal vectors
	<A> foreach 1..($n+5);

	$_=<A>;die "Incorrect header\n" unless /^\t\t(\d*),$/;
	$q=$1;

	# Read in the triangle face information
	foreach $i (1..$q) {
		$_=<A>;
		/^\t\t<(\d*),(\d*),(\d*)>,?$/;
		$b1=$1+1+$vs[$l];
		$b2=$2+1+$vs[$l];
		$b3=$3+1+$vs[$l];
		#print "$b1 $b2 $b3\n";
		vert($b1,$b2);
		vert($b2,$b3);
		vert($b3,$b1);
	}

	$l++;
	$vs[$l]=$vs[$l-1]+$n;
	<A>;<A>;<A>;
}
close A;

# Save the output file
open B,">$ARGV[1]";
if($opt_g) {
	foreach $ll (1..$l) {
		print B "# [Cell $ll]\n";
		print "# [Cell $ll]\n";
		foreach $i (($vs[$ll-1]+1)..$vs[$ll]) {
			$l=1;
			print "$i $v[$i][0]\n";
			while($l+1<=$v[$i][0]) {
				print "here\n";
				$q=$v[$i][$l++];
				$q2=$v[$i][$l++];
				print B "$x[$q] $y[$q] $z[$q]\n$x[$i] $y[$i] $z[$i]\n$x[$q2] $y[$q2] $z[$q2]\n\n";
			}
			if($l==$v[$i][0]) {
				$q=$v[$i][$l];
				print B "$x[$q] $y[$q] $z[$q]\n$x[$i] $y[$i] $z[$i]\n\n";
			}
		}
	}
} else {
	print B "union {\n" if $opt_u;
	foreach $i (1..$n) {
		$q=$v[$i][$_], print B "cylinder{<$x[$i],$y[$i],$z[$i]>,<$x[$q],$y[$q],$z[$q]>,r}\n" foreach 1..$v[$i][0];
	}
	print B "}\n" if $opt_u;
}

sub vert {
	if(@_[0]<@_[1]) {$a=@_[0];$b=@_[1];}
	elsif(@_[0]>@_[1]) {$a=@_[1];$b=@_[0];}
	else {print "Warning: degenerate mesh\n";return;}
	foreach (1..$v[$a][0]) {
		return if $v[$a][$_]==$b;
	}
	#print "-> $a $b\n";
	$v[$a][0]++;$v[$a][$v[$a][0]]=$b;
}

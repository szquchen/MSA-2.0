#!/usr/bin/perl -w
use strict;

die " Syntax:\n   ./postmsa [max degree] [molecule]\n" unless $#ARGV > 0;

my $degree = shift(@ARGV);
my @atom = @ARGV;

my $stub = "MOL_".join("_",@atom)."_$degree";
my $mono = $stub.".MONO";
my $poly = $stub.".POLY";
my $fbas = "basis.f90";


# ... read in the MONO data file
my @datMono = ();
open(MONO,$mono)  or die "Can not open $mono to read!\n";
while(<MONO>){
    chomp;
    my @mono = split;
    push(@datMono,\@mono);
}
close(MONO);

my $nx = $#{$datMono[0]};
my $nm = $#datMono;

# ... read in the POLY data file
my @datPoly = ();
open(POLY,$poly)  or die "Can not open $poly to read!\n";
while(<POLY>){
    chomp;
    my @poly = split;
    push(@datPoly,\@poly);
}
close(POLY);

my $np = $#datPoly;

# ... record the basis function bemsa

open(OUT,">$fbas") or die "Can not open $fbas to write!\n";
print_module_head(\*OUT,$nx,$nm,$np);


print_mono_head(\*OUT,$nx,$nm);
for(my $i=0;$i<=$nm;$i++){print_mono_body(\*OUT, $i, $datMono[$i]);}
print_mono_foot(\*OUT);

print_poly_head(\*OUT, $nm, $np);
for(my $i=0;$i<=$np;$i++){print_poly_body(\*OUT, $i, $datPoly[$i]);}
print_poly_foot(\*OUT);
print_module_foot(\*OUT);
close(OUT);

# ... Here are some subroutines

sub print_module_head{
    my ($out,$nx,$nm,$np) = @_;
    print $out "module basis\n";
    print $out "  implicit none\n";
    print $out "\n";
    print $out "contains\n";
    print $out "  function emsav(x,c) result(v)\n";
    print $out "    implicit none\n";
    print $out "    real,dimension(1:$nx)::x\n";
    print $out "    real,dimension(0:$np)::c\n";
    print $out "    real::v\n";
    print $out "    ! ::::::::::::::::::::\n";
    print $out "    real,dimension(0:$np)::p\n";
    print $out "\n";
    print $out "    call bemsav(x,p)\n";
    print $out "    v = dot_product(p,c)\n";
    print $out "\n";
    print $out "    return\n";
    print $out "  end function emsav\n";
    print $out "\n";
    print $out "  subroutine bemsav(x,p)\n";
    print $out "    implicit none\n";
    print $out "    real,dimension(1:$nx),intent(in)::x\n";
    print $out "    real,dimension(0:$np),intent(out)::p\n";
    print $out "    ! ::::::::::::::::::::\n";
    print $out "    real,dimension(0:$nm)::m\n";
    print $out "\n";
    print $out "    call evmono(x,m)\n";
    print $out "    call evpoly(m,p)\n";
    print $out "\n";
    print $out "    return\n";
    print $out "  end subroutine bemsav\n";
    print $out "\n";
}

sub print_module_foot{
    my ($out) = @_;
    print $out "end module basis\n";
}


sub print_mono_head{
    my ($om, $nx, $nm)= @_;

    print $om "  subroutine evmono(x,m)\n";
    print $om "    implicit none\n";
    print $om "    real,dimension(1:$nx),intent(in)::x\n";
    print $om "    real,dimension(0:$nm),intent(out)::m\n";
    print $om "    !::::::::::::::::::::\n";
    print $om "\n";
    
}

sub print_poly_head{
    my ($op, $nm, $np)= @_;

    print $op "  subroutine evpoly(m,p)\n";
    print $op "    implicit none\n";
    print $op "    real,dimension(0:$nm),intent(in)::m\n";
    print $op "    real,dimension(0:$np),intent(out)::p\n";
    print $op "    !::::::::::::::::::::\n";
    print $op "\n";
    
}

sub print_mono_foot{
    my ($om)= @_;
    
    print $om "\n";
    print $om "    return\n";
    print $om "  end subroutine evmono\n";
    print $om "\n";
    
}

sub print_poly_foot{
    my ($op)= @_;
    
    print $op "\n";
    print $op "    return\n";
    print $op "  end subroutine evpoly\n";
    print $op "\n";
}

sub print_mono_body{
    my ($om, $iMono, $mono)= @_;
    
    printf $om "    m(%d) = ", $iMono;
    
    my $stat = shift(@$mono);
    if($stat == 0){
	my $new = 1;
	my $nz = 0;
	for(my $i=0;$i<=$#$mono;$i++){
	    my $deg = $mono->[$i];
	    next if $deg == 0;
	    $nz++;
	    if($deg == 1){
		$new? printf $om "x(%d)",$i+1: printf $om "*x(%d)",$i+1;
		$new = 0;
	    }else{
		$new? printf $om "x(%d)**$deg",$i+1 : printf $om "*x(%d)**$deg",$i+1;
		$new = 0;
	    }
	    print $om " &\n\t& " if $nz % 5 == 0 && $i<$#$mono;
	}
	print $om "1.0D0" if $new;
    }else{
	printf $om "m(%d)",$mono->[0];
	for(my $i=1;$i<=$#$mono;$i++){
	    my $idx = $mono->[$i];
	    printf $om "*m(%d)",$idx;
	    print $om " &\n\t& " if $i % 5 == 0 && $i<$#$mono;
	}
    }
    print $om "\n";
    unshift(@$mono,$stat);
}

sub print_poly_body{
    my ($op, $iPoly, $poly)= @_;
    
    my $stat = shift(@$poly);
    printf $op "    p(%d) = ",$iPoly;
    if($stat == 2){		# sum of mono
	printf $op "m(%d)",$poly->[0];
	for(my $i=1;$i<=$#$poly;$i++){
	    my $idx = $poly->[$i];
	    printf $op " + m(%d)",$idx;
	    print $op " &\n\t& " if $i % 5 == 0 && $i<$#$poly;
	}
    }else{			# usable decomp
	my $ima = shift(@$poly);
	my $imb = shift(@$poly);
	printf $op "p(%d)*p(%d)", $ima, $imb;
	for(my $i=0;$i<=$#$poly;$i++){
	    my $idx = $poly->[$i];
	    printf $op " - p(%d)",$idx;
	    print $op " &\n\t& " if $i % 5 == 0 && $i<$#$poly && $i > 0;
	}
	
	unshift(@$poly,$imb);
	unshift(@$poly,$ima);
    }
    print $op "\n";
    unshift(@$poly,$stat);
}

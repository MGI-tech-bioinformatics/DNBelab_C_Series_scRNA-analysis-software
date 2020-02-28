#!/usr/bin/perl
use strict;
use warnings;
use List::Util;
use Statistics::Descriptive;
use Getopt::Long;
my($input, $output, $length, $help, $threshold, $min, $num, $limit);
my%hash=(
    0 => 'A',
    1 => 'T',
    2 => 'C',
    3 => 'G',
);
GetOptions(
    'i|input:s' => \$input,
    'o|output:s' => \$output,
    'l|length:i' => \$length,
    'h|help:s' => \$help,
    't|threshold:i' => \$threshold,
    'm|min:i' => \$min,
    'n|num:i' => \$num,
    'limit:i' => \$limit,
);
sub useage {
	print"\nusage $0 -i <input> -o <output> -l <length>";
	my $text_main="
	 -i	Input file
	 -o	Output file
	 -l	Barcode length, default[20]
	 -t	Threshold of colorbalance, default[0.125]
	 -m	Minimum barcodes number, default[96]
         -n     Shuffle number, default[5]
         -limit Cycle times, default[10]
	 -h	Help
	 ";
	 print $text_main,"\n";
	 exit();
 }
if(!$input||!$output){
    useage();
}
$length ||= 20;
$threshold ||= 0.125;
$min ||= 96;
$num ||= 5;
$limit ||=10;
 
 #Calculate the base composition in each position
#计算输入的序列集中，指定长度范围0-len-1的碱基比例 
 sub BiasCount {
    my $input_ref =shift;
    my $len = shift;
    my @input = @$input_ref;
    my (%A, %T, %C, %G);
    my $ref = [\%A, \%T, \%C, \%G];
    my $line = @input;
    my %tmp;
    for my $i (0...($len-1)){
        ($A{$i}, $T{$i}, $C{$i}, $G{$i}) = (0, 0, 0, 0);
        foreach (@input){
            chomp;
            my @cache = split;
            $tmp{$cache[0]} = $cache[1];
            my $seq = $cache[1];
            my $tmp = substr($seq, $i, 1);
            if ($tmp =~ tr/A//){
                $A{$i} += 1;
            }
            elsif($tmp =~ tr/T//){
                $T{$i} += 1;
            }
            elsif($tmp =~ tr/C//){
                $C{$i} += 1;
            }
            elsif($tmp =~ tr/G//){
                $G{$i} += 1;
            }
        }
    }
    foreach $a (0...($len-1)){
        $A{$a} = ($A{$a} / $line);
        $T{$a} = ($T{$a} / $line);
        $C{$a} = ($C{$a} / $line);
        $G{$a} = ($G{$a} / $line);
    }
    $input_ref = 0;
    return $ref, $line, %tmp;
};

#Get the min bias index  in each position

sub FindIndex {
    my $ref = shift;
    my $a = shift;
    my $bias;
    my @bias = ($ref->[0]{$a}, $ref->[1]{$a}, $ref->[2]{$a}, $ref->[3]{$a});
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(\@bias);
    my $mindex = $stat->mindex();
    $ref = 0;
    return $bias = $hash{$mindex};
}

#Counting frequencies by different lengths of barcodes

sub MinLengthIndex {
    my $input_ref = shift;
    my @input = @$input_ref;
    my%hash;
    my@len;
    my@min_idx;
    foreach (@input){
        chomp;
	    #    [index]	[seq]
	my @tmp = split;
	my $seq = (split /AATTC/, $tmp[1])[0]; 
	$hash{length $seq} += 1;
    };
    foreach my$key(sort{$hash{$a}<=>$hash{$b}}keys%hash){
        push @len, $key;
    };
    foreach(@input){
	chomp;
	my @tmp = split;
	my $seq = (split /AATTC/, $tmp[1])[0]; 
        foreach my$len(@len[0...3]){
	    if(length($seq) == $len){
	        push @min_idx, $tmp[0];
	    };
	};
    };
    $input_ref = 0;
    return @min_idx;
};

#Counting the false index

sub MakeIndexList {
    my $input_ref = shift;
    my $pos_ref = shift;
    my $bias_ref = shift;
    my @input = @$input_ref;
    my @pos = @$pos_ref;
    my @bias = @$bias_ref;
    my @c;
    my @out;
    foreach (@input) {
	chomp;
		#    [index]	[seq]
	my @tmp = split;
	for my $i (0...$#pos) {
	    my $s = substr($tmp[1], $pos[$i], 1);
	    my $c = ($s eq $bias[$i])?1:0;
	    push @c, $c;
	};
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data(\@c);
        my $sum = $stat->sum();
        if ($sum == 0){
            push @out, $tmp[0];
        };
        @c = ();
    };
    ($input_ref, $pos_ref, $bias_ref) = (0, 0, 0);
    return @out;
};

sub RandemFirstLength {
    my $input_ref = shift;
    my $tmp_ref = shift;
    my @input = @$input_ref;
    my %tmp = %$tmp_ref;
    my%hash;
    my(@len, @seq);
    my$firstlen;
    foreach (@input){
        chomp;
        my @tmp = split;
        my $seq = (split /AATTC/, $tmp[1])[0];
        $hash{length $seq} += 1;
    };
    foreach my$key(sort{$hash{$b}<=>$hash{$a}}keys%hash){
        push @len, $key;
    };
    $firstlen = shift @len;
    my @keys = keys %tmp;
    @keys = List::Util::shuffle @keys;
    for my$i(0...$#keys){
        my $ss = (split /AATTC/, $tmp{$keys[$i]})[0];
        if (length $ss == $firstlen){
            splice @keys, $i, 1;
            last;
        }
    }
    @input = ();
    foreach (@keys){
        push @input, ($_."\t".$tmp{$_});
    }
    ($input_ref, $tmp_ref) = (0, 0);
    return @input;
}


################################################################
#                                                              #
#                 Start to do some magic !                     #
#                                                              #
################################################################

my $time = 1;
my@in;
my@input;
my$false = 0;
open INPUT, "$input" or die print "Can't open $input","\n" && useage();
foreach (<INPUT>){
	chomp;
	push @in, $_;
};

print "Star the colorbalance!!\n";
close INPUT;
while ($time < ($limit+1)){
@input = @in;
while (1){
    my @bias = ();
    my @pos = ();
    my %tmp;
    my ($ref, $line);
    ($ref, $line, %tmp) = BiasCount(\@input, $length);   
    for my $a (0...($length-1)){
        if ($ref->[0]{$a} < $threshold || $ref->[1]{$a} < $threshold || $ref->[2]{$a} < $threshold || $ref->[3]{$a} < $threshold){
	    my $bias = FindIndex($ref, $a);
	    push @bias, $bias;
	    push @pos, $a;
	}
    };
    if(@bias == 0 && $line > $min){
        @input = RandemFirstLength(\@input, \%tmp);
        if(15 >= ($line - $min)){
            $num = 1;
        }
        redo;
    }
    elsif(@bias == 0 && $line == $min){
        $false = 0;
        last;
    }
    if ($min > $line||@bias == ($length/2)){
        $false = 1;
        last;
    } 
    my @out = MakeIndexList(\@input, \@pos, \@bias);   # Get the filter out barcode index
    my @min_idx = MinLengthIndex(\@input);           # Get the save barcode index
    my %out = map{$_ => 1}@out;
    my %min_idx = map{$_ => 1}@min_idx;
    my @out_only = grep{!$min_idx{$_}}@out;
    @out_only=List::Util::shuffle @out_only;
    if (@out_only == 0){
        $false = 1;
        last;
    }
    my %out_only = map{$_ => 1}@out_only[0...($num-1)];
    my @key = sort keys %tmp;
    my @key_filter = grep{!$out_only{$_}}@key;
    @input=();
    foreach (@key_filter){
        push @input, ($_."\t".$tmp{$_});
    }
}
print "Let's starting $time times !\n";
$time += 1;
if ($time == ($limit+1)){
    print "\nSorry! Barcode list can't be colorbalance, please add limit number or barcodelist!\n";
}
if ($false == 0){
    last;
}
}


my $index = 1;
open OUTPUT, ">$output" or die print "Can't open $output","\n" && useage();  
foreach (@input){
    chomp;
    my $seq = (split)[1];
    print OUTPUT $index."\t".$seq."\n";
    $index += 1;
}
my $total = $index - 1;
print "\nLong live the BGI ! We are success, $total barcodes have been colorbalance.\n" if ($false == 0);
close OUTPUT;	
	
#############################################################################
#                                                                           # 
#                          Long live the BGI                                #
#                                                                           #
#############################################################################

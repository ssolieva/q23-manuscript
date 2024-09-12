#!/usr/bin/perl

# OAS does not add donor or search critera to every sequence line, which we need to do if we want to split sequences evenly over files and recover information such as donor.

my %fields; # additional field container extracted from header line of each file.
open(TMP, "$ARGV[0]") or die "No such file $ARGV[0]";
while ($line = <TMP>){
    $line =~ s/\n//;
    if ($line =~ /Run/){
#	print "HEADER\n";

	# Remove brackets
	$line =~ s/[{}]//g;
	#	print $line."\n";
	
	# Add double quotes to numbered fields
	$line =~ s/\"\"(.+?)\"\"\: ([0-9]+)/\"\"$1\"\": \"\"$2\"\"/g;
	#	print $line."\n";

	# Subsitute field delimiter ',' for '+-*' (because people use commas within strings so hard to parse using it)
	$line =~ s/\"\"(.+?)\"\"\: \"\"(.+?)\"\",/$1:$2+-*/g;

	# Get rid of extra quotes
	$line =~ s/\"//g;
	#	print $line."\n";

	# Tokenize by new field delimiter
	my @toks = split /\+\-\*/,$line;

	#  Split into key->value pairs
	foreach (@toks){
	    my @pairs = split(":",$_);
	    # Remove space at start
	    $pairs[0] =~ s/^ (.+)/$1/;
	    $pairs[1] =~ s/^ (.+)/$1/;

	    # Remove stupid commas
	    $pairs[0] =~ s/\,//g;
	    $pairs[1] =~ s/\,//g;
	    
	    #print "$pairs[0] = $pairs[1]\n";
	    $fields{$pairs[0]} = $pairs[1];
	    
	}
	next;
    }

    # Normal field description line
    if ($line =~ /^sequence/){
	$line = "Run,Subject,Author,Vaccine,Age,BType,Disease,Isotype,".$line;
	print $line."\n";
	next;
    }

    # A data line
    if ($line =~ /^[AGTC]/){
	$line = $fields{"Run"}.",".$fields{"Subject"}.",".$fields{"Author"}.",".$fields{"Vaccine"}.",".$fields{"Age"}.",".$fields{"BType"}.",".$fields{"Disease"}.",".$fields{"Isotype"}.",".$line;
	print $line."\n";
    }
}

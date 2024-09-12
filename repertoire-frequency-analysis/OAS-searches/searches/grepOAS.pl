#!/bin/perl
use File::Basename;
use Getopt::Long;
use String::LCSS_XS qw(lcss lcss_all);

GetOptions ("file=s" => \$file,
	    "format=s" =>\$format,
	    "keepX" =>\$keepX,  # Keep sequences with Xs in them
	    
	    # Regex for whole HC or LC
	    "regexHC=s" => \$regexHC,
	    "regexLC=s" => \$regexLC,
	    
	    "regexHCDR1=s" => \$regexHCDR1,            
	    "regexHCDR2=s" => \$regexHCDR2,            
	    "regexHCDR3=s" => \$regexHCDR3,            

	    "regexLCDR1=s" => \$regexLCDR1,            
	    "regexLCDR2=s" => \$regexLCDR2,            
	    "regexLCDR3=s" => \$regexLCDR3,            
	    
	    "vhgene=s" => \$VHgene,
	    "dhgene=s" => \$DHgene,
	    "jhgene=s" => \$JHgene,

	    "vlgene=s" => \$VLgene,
	    "jlgene=s" => \$JLgene,
	    
            "Hcdr1len=s" =>\$Hcdr1len,
	    "Hcdr2len=s" =>\$Hcdr2len,
	    "Hcdr3len=s" =>\$Hcdr3len,
	    "Lcdr1len=s" =>\$Lcdr1len,
	    "Lcdr2len=s" =>\$Lcdr2len,
	    "Lcdr3len=s" =>\$Lcdr3len,
	    
	    "outFasta=s" =>\$outFasta,  # Output a fasta file with named with this string
	    "seqFasta=s" =>\$fasta_seq); # use field name to populate fasta

# Default values
$jobname="grepCDR3" if (!$jobname);
$regexHC = 1 if (!$regexHC);
$regexLC = 1 if (!$regexLC);
$regexHCDR1 = 1 if (!$regexHCDR1);
$regexHCDR2 = 1 if (!$regexHCDR2);
$regexHCDR3 = 1 if (!$regexHCDR3);

$regexLCDR1 = 1 if (!$regexLCDR1);
$regexLCDR2 = 1 if (!$regexLCDR2);
$regexLCDR3 = 1 if (!$regexLCDR3);

$VHgene = 1 if (!$VHgene);
$DHgene = 1 if (!$DHgene);
$JHgene = 1 if (!$JHgene);
$VLgene = 1 if (!$VLgene);
$JLgene = 1 if (!$JLgene);

$format = "OASheavy" if (!$format);

my %codons_all = (
'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
'GTG' => 'V',
'TAA' => 'STOP',
'TGA' => 'STOP',
'TAG' => 'STOP',);

my $Hcdr1len_small = 0;
my $Hcdr1len_big   = 100;
if (!$Hcdr1len) {
    $Hcdr1len = 1;
} else {
    if ($Hcdr1len !~ /:/) { die ("CDRLen must be in format START:END,e.g 10:25\n"); }
    my @toks  = split(":",$Hcdr1len);
    $Hcdr1len_small = $toks[0];
    $Hcdr1len_big   = $toks[1];
}
my $Hcdr2len_small = 0;
my $Hcdr2len_big   = 100;
if (!$Hcdr2len) {
    $Hcdr2len = 1;
} else {
    if ($Hcdr2len !~ /:/) { die ("CDRLen must be in format START:END,e.g 10:25\n"); }
    my @toks  = split(":",$Hcdr2len);
    $Hcdr2len_small = $toks[0];
    $Hcdr2len_big   = $toks[1];
}
my $Hcdr3len_small = 0;
my $Hcdr3len_big   = 100;
if (!$Hcdr3len) {
    $Hcdr3len = 1;
} else {
    if ($Hcdr3len !~ /:/) { die ("CDRLen must be in format START:END,e.g 10:25\n"); }
    my @toks  = split(":",$Hcdr3len);
    $Hcdr3len_small = $toks[0];
    $Hcdr3len_big   = $toks[1];
}

my $Lcdr1len_small = 0;
my $Lcdr1len_big   = 100;
if (!$Lcdr1len) {
    $Lcdr1len = 1;
} else {
    if ($Lcdr1len !~ /:/) { die ("CDRLen must be in format START:END,e.g 10:25\n"); }
    my @toks  = split(":",$Lcdr1len);
    $Lcdr1len_small = $toks[0];
    $Lcdr1len_big   = $toks[1];
}
my $Lcdr2len_small = 0;
my $Lcdr2len_big   = 100;
if (!$Lcdr2len) {
    $Lcdr2len = 1;
} else {
    if ($Lcdr2len !~ /:/) { die ("CDRLen must be in format START:END,e.g 10:25\n"); }
    my @toks  = split(":",$Lcdr2len);
    $Lcdr2len_small = $toks[0];
    $Lcdr2len_big   = $toks[1];
}
my $Lcdr3len_small = 0;
my $Lcdr3len_big   = 100;
if (!$Lcdr3len) {
    $Lcdr3len = 1;
} else {
    if ($Lcdr3len !~ /:/) { die ("CDRLen must be in format START:END,e.g 10:25\n"); }
    my @toks  = split(":",$Lcdr3len);
    $Lcdr3len_small = $toks[0];
    $Lcdr3len_big   = $toks[1];
}



my %names;
$names{"seq_id"}=0;

my $header = "";
# modifed the following to match abstar output
# seq_id -> there isn't any  
# v_full -> v_call  (v_full_heavy)
# d_full -> d_call  (d_full_heavy)
# j_full -> j_call  (j_full_heavy)
# cdr3_length -> junction_aa_length (cdr3_length)
# cdr3_nt -> cdr3
# cdr3_aa -> cdr3_aa
# cdr3_length = junction_aa_length_heavy
# v_full      = v_call_heavy
# d_full      = d_call_heavy
# j_full      = j_call_heavy

if ($format eq "OASheavy"){
    #    $header="Run,Subject,Author,Vaccine,Age,BType,Disease,Isotype,sequence,locus,stop_codon,vj_in_frame,v_frameshift,productive,rev_comp,complete_vdj,v_call,d_call,j_call,sequence_alignment,germline_alignment,sequence_alignment_aa,germline_alignment_aa,v_alignment_start,v_alignment_end,d_alignment_start,d_alignment_end,j_alignment_start,j_alignment_end,v_sequence_alignment,v_sequence_alignment_aa,v_germline_alignment,v_germline_alignment_aa,d_sequence_alignment,d_sequence_alignment_aa,d_germline_alignment,d_germline_alignment_aa,j_sequence_alignment,j_sequence_alignment_aa,j_germline_alignment,j_germline_alignment_aa,fwr1,fwr1_aa,cdr1,cdr1_aa,fwr2,fwr2_aa,cdr2,cdr2_aa,fwr3,fwr3_aa,fwr4,fwr4_aa,cdr3,cdr3_aa,junction,junction_length,junction_aa,junction_aa_length,v_score,d_score,j_score,v_cigar,d_cigar,j_cigar,v_support,d_support,j_support,v_identity,d_identity,j_identity,v_sequence_start,v_sequence_end,v_germline_start,v_germline_end,d_sequence_start,d_sequence_end,d_germline_start,d_germline_end,j_sequence_start,j_sequence_end,j_germline_start,j_germline_end,fwr1_start,fwr1_end,cdr1_start,cdr1_end,fwr2_start,fwr2_end,cdr2_start,cdr2_end,fwr3_start,fwr3_end,fwr4_start,fwr4_end,cdr3_start,cdr3_end,np1,np1_length,np2,np2_length,c_region,Redundancy,ANARCI_numbering,ANARCI_status";
    $header="Run,Subject,Author,Vaccine,Age,BType,Disease,Isotype,sequence_heavy,locus_heavy,stop_codon_heavy,vj_in_frame_heavy,v_frameshift_heavy,productive_heavy,rev_comp_heavy,complete_vdj_heavy,v_call_heavy,d_call_heavy,j_call_heavy,sequence_alignment_heavy,germline_alignment_heavy,sequence_alignment_aa_heavy,germline_alignment_aa_heavy,v_alignment_start_heavy,v_alignment_end_heavy,d_alignment_start_heavy,d_alignment_end_heavy,j_alignment_start_heavy,j_alignment_end_heavy,v_sequence_alignment_heavy,v_sequence_alignment_aa_heavy,v_germline_alignment_heavy,v_germline_alignment_aa_heavy,d_sequence_alignment_heavy,d_sequence_alignment_aa_heavy,d_germline_alignment_heavy,d_germline_alignment_aa_heavy,j_sequence_alignment_heavy,j_sequence_alignment_aa_heavy,j_germline_alignment_heavy,j_germline_alignment_aa_heavy,fwr1_heavy,fwr1_aa_heavy,cdr1_heavy,cdr1_aa_heavy,fwr2_heavy,fwr2_aa_heavy,cdr2_heavy,cdr2_aa_heavy,fwr3_heavy,fwr3_aa_heavy,fwr4_heavy,fwr4_aa_heavy,cdr3_heavy,cdr3_aa_heavy,junction_heavy,junction_length_heavy,junction_aa_heavy,junction_aa_length_heavy,v_score_heavy,d_score_heavy,j_score_heavy,v_cigar_heavy,d_cigar_heavy,j_cigar_heavy,v_support_heavy,d_support_heavy,j_support_heavy,v_identity_heavy,d_identity_heavy,j_identity_heavy,v_sequence_start_heavy,v_sequence_end_heavy,v_germline_start_heavy,v_germline_end_heavy,d_sequence_start_heavy,d_sequence_end_heavy,d_germline_start_heavy,d_germline_end_heavy,j_sequence_start_heavy,j_sequence_end_heavy,j_germline_start_heavy,j_germline_end_heavy,fwr1_start_heavy,fwr1_end_heavy,cdr1_start_heavy,cdr1_end_heavy,fwr2_start_heavy,fwr2_end_heavy,cdr2_start_heavy,cdr2_end_heavy,fwr3_start_heavy,fwr3_end_heavy,fwr4_start_heavy,fwr4_end_heavy,cdr3_start_heavy,cdr3_end_heavy,np1_heavy,np1_length_heavy,np2_heavy,np2_length_heavy,c_region_heavy,Redundancy_heavy,ANARCI_numbering_heavy,ANARCI_status_heavy";
}

if ($format eq "OASlight"){
    #    $header="Run,Subject,Author,Vaccine,Age,BType,Disease,Isotype,sequence,locus,stop_codon,vj_in_frame,v_frameshift,productive,rev_comp,complete_vdj,v_call,d_call,j_call,sequence_alignment,germline_alignment,sequence_alignment_aa,germline_alignment_aa,v_alignment_start,v_alignment_end,d_alignment_start,d_alignment_end,j_alignment_start,j_alignment_end,v_sequence_alignment,v_sequence_alignment_aa,v_germline_alignment,v_germline_alignment_aa,d_sequence_alignment,d_sequence_alignment_aa,d_germline_alignment,d_germline_alignment_aa,j_sequence_alignment,j_sequence_alignment_aa,j_germline_alignment,j_germline_alignment_aa,fwr1,fwr1_aa,cdr1,cdr1_aa,fwr2,fwr2_aa,cdr2,cdr2_aa,fwr3,fwr3_aa,fwr4,fwr4_aa,cdr3,cdr3_aa,junction,junction_length,junction_aa,junction_aa_length,v_score,d_score,j_score,v_cigar,d_cigar,j_cigar,v_support,d_support,j_support,v_identity,d_identity,j_identity,v_sequence_start,v_sequence_end,v_germline_start,v_germline_end,d_sequence_start,d_sequence_end,d_germline_start,d_germline_end,j_sequence_start,j_sequence_end,j_germline_start,j_germline_end,fwr1_start,fwr1_end,cdr1_start,cdr1_end,fwr2_start,fwr2_end,cdr2_start,cdr2_end,fwr3_start,fwr3_end,fwr4_start,fwr4_end,cdr3_start,cdr3_end,np1,np1_length,np2,np2_length,c_region,Redundancy,ANARCI_numbering,ANARCI_status";
    $header="Run,Subject,Author,Vaccine,Age,BType,Disease,Isotype,sequence_light,locus_light,stop_codon_light,vj_in_frame_light,v_frameshift_light,productive_light,rev_comp_light,complete_vdj_light,v_call_light,d_call_light,j_call_light,sequence_alignment_light,germline_alignment_light,sequence_alignment_aa_light,germline_alignment_aa_light,v_alignment_start_light,v_alignment_end_light,d_alignment_start_light,d_alignment_end_light,j_alignment_start_light,j_alignment_end_light,v_sequence_alignment_light,v_sequence_alignment_aa_light,v_germline_alignment_light,v_germline_alignment_aa_light,d_sequence_alignment_light,d_sequence_alignment_aa_light,d_germline_alignment_light,d_germline_alignment_aa_light,j_sequence_alignment_light,j_sequence_alignment_aa_light,j_germline_alignment_light,j_germline_alignment_aa_light,fwr1_light,fwr1_aa_light,cdr1_light,cdr1_aa_light,fwr2_light,fwr2_aa_light,cdr2_light,cdr2_aa_light,fwr3_light,fwr3_aa_light,fwr4_light,fwr4_aa_light,cdr3_light,cdr3_aa_light,junction_light,junction_length_light,junction_aa_light,junction_aa_length_light,v_score_light,d_score_light,j_score_light,v_cigar_light,d_cigar_light,j_cigar_light,v_support_light,d_support_light,j_support_light,v_identity_light,d_identity_light,j_identity_light,v_sequence_start_light,v_sequence_end_light,v_germline_start_light,v_germline_end_light,d_sequence_start_light,d_sequence_end_light,d_germline_start_light,d_germline_end_light,j_sequence_start_light,j_sequence_end_light,j_germline_start_light,j_germline_end_light,fwr1_start_light,fwr1_end_light,cdr1_start_light,cdr1_end_light,fwr2_start_light,fwr2_end_light,cdr2_start_light,cdr2_end_light,fwr3_start_light,fwr3_end_light,fwr4_start_light,fwr4_end_light,cdr3_start_light,cdr3_end_light,np1_light,np1_length_light,np2_light,np2_length_light,c_region_light,Redundancy_light,ANARCI_numbering_light,ANARCI_status_light";
}

if ($format eq "OASpaired"){
    
    $header = "Run,Subject,Author,Vaccine,Age,BType,Disease,Isotype,sequence_id_heavy,sequence_heavy,locus_heavy,stop_codon_heavy,vj_in_frame_heavy,productive_heavy,rev_comp_heavy,v_call_heavy,d_call_heavy,j_call_heavy,sequence_alignment_heavy,germline_alignment_heavy,sequence_alignment_aa_heavy,germline_alignment_aa_heavy,v_alignment_start_heavy,v_alignment_end_heavy,d_alignment_start_heavy,d_alignment_end_heavy,j_alignment_start_heavy,j_alignment_end_heavy,v_sequence_alignment_heavy,v_sequence_alignment_aa_heavy,v_germline_alignment_heavy,v_germline_alignment_aa_heavy,d_sequence_alignment_heavy,d_sequence_alignment_aa_heavy,d_germline_alignment_heavy,d_germline_alignment_aa_heavy,j_sequence_alignment_heavy,j_sequence_alignment_aa_heavy,j_germline_alignment_heavy,j_germline_alignment_aa_heavy,fwr1_heavy,fwr1_aa_heavy,cdr1_heavy,cdr1_aa_heavy,fwr2_heavy,fwr2_aa_heavy,cdr2_heavy,cdr2_aa_heavy,fwr3_heavy,fwr3_aa_heavy,cdr3_heavy,cdr3_aa_heavy,junction_heavy,junction_length_heavy,junction_aa_heavy,junction_aa_length_heavy,v_score_heavy,d_score_heavy,j_score_heavy,v_cigar_heavy,d_cigar_heavy,j_cigar_heavy,v_support_heavy,d_support_heavy,j_support_heavy,v_identity_heavy,d_identity_heavy,j_identity_heavy,v_sequence_start_heavy,v_sequence_end_heavy,v_germline_start_heavy,v_germline_end_heavy,d_sequence_start_heavy,d_sequence_end_heavy,d_germline_start_heavy,d_germline_end_heavy,j_sequence_start_heavy,j_sequence_end_heavy,j_germline_start_heavy,j_germline_end_heavy,fwr1_start_heavy,fwr1_end_heavy,cdr1_start_heavy,cdr1_end_heavy,fwr2_start_heavy,fwr2_end_heavy,cdr2_start_heavy,cdr2_end_heavy,fwr3_start_heavy,fwr3_end_heavy,cdr3_start_heavy,cdr3_end_heavy,np1_heavy,np1_length_heavy,np2_heavy,np2_length_heavy,sequence_id_light,sequence_light,locus_light,stop_codon_light,vj_in_frame_light,productive_light,rev_comp_light,v_call_light,d_call_light,j_call_light,sequence_alignment_light,germline_alignment_light,sequence_alignment_aa_light,germline_alignment_aa_light,v_alignment_start_light,v_alignment_end_light,d_alignment_start_light,d_alignment_end_light,j_alignment_start_light,j_alignment_end_light,v_sequence_alignment_light,v_sequence_alignment_aa_light,v_germline_alignment_light,v_germline_alignment_aa_light,d_sequence_alignment_light,d_sequence_alignment_aa_light,d_germline_alignment_light,d_germline_alignment_aa_light,j_sequence_alignment_light,j_sequence_alignment_aa_light,j_germline_alignment_light,j_germline_alignment_aa_light,fwr1_light,fwr1_aa_light,cdr1_light,cdr1_aa_light,fwr2_light,fwr2_aa_light,cdr2_light,cdr2_aa_light,fwr3_light,fwr3_aa_light,cdr3_light,cdr3_aa_light,junction_light,junction_length_light,junction_aa_light,junction_aa_length_light,v_score_light,d_score_light,j_score_light,v_cigar_light,d_cigar_light,j_cigar_light,v_support_light,d_support_light,j_support_light,v_identity_light,d_identity_light,j_identity_light,v_sequence_start_light,v_sequence_end_light,v_germline_start_light,v_germline_end_light,d_sequence_start_light,d_sequence_end_light,d_germline_start_light,d_germline_end_light,j_sequence_start_light,j_sequence_end_light,j_germline_start_light,j_germline_end_light,fwr1_start_light,fwr1_end_light,cdr1_start_light,cdr1_end_light,fwr2_start_light,fwr2_end_light,cdr2_start_light,cdr2_end_light,fwr3_start_light,fwr3_end_light,cdr3_start_light,cdr3_end_light,np1_light,np1_length_light,np2_light,np2_length_light,ANARCI_numbering_light,ANARCI_numbering_heavy,ANARCI_status_light,ANARCI_status_heavy";
}

my @fields = split(",",$header);
for ($i=0; $i<$#fields;$i++){
    $names{$fields[$i]} = $i;
}

die "Need a file to run grepOAS.pl\n" if (!$file);

# Basename of file
my $bfile = fileparse($file);

if ($outFasta){
    open(FAS, '>',$bfile."_".$outFasta);
}

my $line_count = 0;
open(TMP,"<","$file") or die "Couldn't open $file\n";
while ($line = <TMP>){

    # Remove whitespace
   $line =~ s/\n//;

   # Tokenize line on comma
   my @toks = split(",",$line);

   # Default to removing sequences with 'X's in them
   next if (!$keepX && $format eq "OASheavy" && $toks[$names{"sequence_alignment_aa_heavy"}] =~ /X/);
   next if (!$keepX && $format eq "OASlight" && $toks[$names{"sequence_alignment_aa_light"}] =~ /X/);
   next if (!$keepX && $format eq "OASpaired" && ($toks[$names{"sequence_alignment_aa_heavy"}] =~ /X/ || $toks[$names{"sequence_alignment_aa_light"}] =~ /X/));
   
   # Need to modify a few fields based on differences in different databases
   # Store sequence id
    my $seq_id = $toks[$names{"seq_id"}];

   # Store CDR3len
   my $Hcdr1len_seq  = length($toks[$names{"cdr1_aa_heavy"}]);
   my $Hcdr2len_seq  = length($toks[$names{"cdr2_aa_heavy"}]);
   my $Hcdr3len_seq  = length($toks[$names{"cdr3_aa_heavy"}]);

   my $Lcdr1len_seq  = length($toks[$names{"cdr1_aa_light"}]);
   my $Lcdr2len_seq  = length($toks[$names{"cdr2_aa_light"}]);
   my $Lcdr3len_seq  = length($toks[$names{"cdr3_aa_light"}]);
       
    # Add line number onto seq_id for OAS format
    $line_count++;
    if ($format eq "OASheavy"){
	$seq_id = "OASH-".$line_count;
   }

   if ($format eq "OASlight"){
	$seq_id = "OASL-".$line_count;
#	$Lcdr3len_seq -= 2;  # atleast for light chains this is true, double check for heavy chains in OAS
    }
   
   if ($format eq "OASpaired"){
       $seq_id = "OASpair-".$toks[$names{"sequence_id_heavy"}]."-".$toks[$names{"sequence_id_light"}];
   }

   if (
        ($regexHC==1 || $toks[$names{"sequence_alignment_aa_heavy"}] =~ /(\S*)$regexHC(\S*)/) &&
        ($regexLC==1 || $toks[$names{"sequence_alignment_aa_light"}] =~ /(\S*)$regexLC(\S*)/) &&

        ($regexHCDR1==1 || $toks[$names{"cdr1_aa_heavy"}] =~ /(\S*)$regexHCDR1(\S*)/) &&
	($regexHCDR2==1 || $toks[$names{"cdr2_aa_heavy"}] =~ /(\S*)$regexHCDR2(\S*)/) &&
	($regexHCDR3==1 || $toks[$names{"cdr3_aa_heavy"}] =~ /(\S*)$regexHCDR3(\S*)/) &&

        ($regexLCDR1==1 || $toks[$names{"cdr1_aa_light"}] =~ /(\S*)$regexLCDR1(\S*)/) &&
	($regexLCDR2==1 || $toks[$names{"cdr2_aa_light"}] =~ /(\S*)$regexLCDR2(\S*)/) &&
	($regexLCDR3==1 || $toks[$names{"cdr3_aa_light"}] =~ /(\S*)$regexLCDR3(\S*)/) &&

	($Hcdr3len == 1  || ($Hcdr3len_seq >= $Hcdr3len_small && $Hcdr3len_seq <= $Hcdr3len_big)) &&
        ($Hcdr2len == 1  || ($Hcdr2len_seq >= $Hcdr2len_small && $Hcdr2len_seq <= $Hcdr2len_big)) &&
        ($Hcdr1len == 1  || ($Hcdr1len_seq >= $Hcdr1len_small && $Hcdr1len_seq <= $Hcdr1len_big)) &&
	($Lcdr3len == 1  || ($Lcdr3len_seq >= $Lcdr3len_small && $Lcdr3len_seq <= $Lcdr3len_big)) &&
        ($Lcdr2len == 1  || ($Lcdr2len_seq >= $Lcdr2len_small && $Lcdr2len_seq <= $Lcdr2len_big)) &&
	($Lcdr1len == 1  || ($Lcdr1len_seq >= $Lcdr1len_small && $Lcdr1len_seq <= $Lcdr1len_big)) &&

	($VHgene == 1   || $toks[$names{"v_call_heavy"}]  =~ /$VHgene/) &&
	($DHgene == 1   || $toks[$names{"d_call_heavy"}]  =~ /$DHgene/) && 
	($JHgene == 1   || $toks[$names{"j_call_heavy"}]  =~ /$JHgene/) &&

	($VLgene == 1   || $toks[$names{"v_call_light"}]  =~ /$VLgene/) &&
	($JLgene == 1   || $toks[$names{"j_call_light"}]  =~ /$JLgene/)){

       if ($format eq "OASheavy"){
	   printf("%-35s %-20s %-25s %-15s %-15s %15s %-15s %-15s %-50s %2d %s %s %s",
		  $seq_id,
		  $bfile,
		  $toks[$names{"Subject"}],
		  $toks[$names{"v_call_heavy"}],
		  $toks[$names{"d_call_heavy"}],
		  $toks[$names{"j_call_heavy"}],
		  $toks[$names{"cdr1_aa_heavy"}],
		  $toks[$names{"cdr2_aa_heavy"}],
		  $toks[$names{"cdr3_aa_heavy"}],
		  $Hcdr3len_seq,
		  $toks[$names{"cdr3_heavy"}],
		  $toks[$names{"sequence_alignment_aa_heavy"}],
		  $toks[$names{"sequence_heavy"}]);
       }

       if ($format eq "OASlight"){
	   printf("%-35s %-20s %-25s %-15s %-15s %-15s %-15s %-15s %2d %s %s %s",
		  $seq_id,
		  $bfile,
		  $toks[$names{"Subject"}],
		  $toks[$names{"v_call_light"}],
		  $toks[$names{"j_call_light"}],
		  $toks[$names{"cdr1_aa_light"}],
		  $toks[$names{"cdr2_aa_light"}],
		  $toks[$names{"cdr3_aa_light"}],
		  $Lcdr3len_seq,
		  $toks[$names{"cdr3_light"}],
		  $toks[$names{"sequence_alignment_aa_light"}],
		  $toks[$names{"sequence_light"}]);
       }
       if ($format eq "OASpaired"){
	     
	   printf("%-35s %-20s %-25s %-15s %-15s %-15s %-15s %15s %-15s %-15s %35s %2d %-15s %-15s %-15s %2d %s %s %s %s",
		  $seq_id,
		  $bfile,
		  $toks[$names{"Subject"}],
		  $toks[$names{"v_call_heavy"}],
		  $toks[$names{"d_call_heavy"}],
		  $toks[$names{"j_call_heavy"}],
		  $toks[$names{"v_call_light"}],
		  $toks[$names{"j_call_light"}],
		  $toks[$names{"cdr1_aa_heavy"}],
		  $toks[$names{"cdr2_aa_heavy"}],
		  $toks[$names{"cdr3_aa_heavy"}],
		  $Hcdr3len_seq,
		  $toks[$names{"cdr1_aa_light"}],
		  $toks[$names{"cdr2_aa_light"}],
		  $toks[$names{"cdr3_aa_light"}],
		  $Lcdr3len_seq,
		  $toks[$names{"sequence_alignment_aa_heavy"}],
		  $toks[$names{"sequence_heavy"}],
	          $toks[$names{"sequence_alignment_aa_light"}],
		  $toks[$names{"sequence_light"}]);
       }

       print "\n";

       # Output a fasta file?
       if ($outFasta){
	   my $vdj = "";
	   if ($format eq "OASheavy" && !$fasta_seq){
	       $vdj =  $toks[$names{"sequence_alignment_aa_heavy"}];
	   }
	   if ($format eq "OASlight"){
	       if (!$fasta_seq){
		   $vdj = $toks[$names{"sequence_alignment_aa_light"}];
	       } 

	       if (($fasta_seq eq "regexLC" && $regexLC==1) || $toks[$names{"sequence_alignment_aa_light"}] =~ /($regexLC)/){
		   $vdj = $1;
	       }
	       
	   }
	   if ($format eq "OASpaired" && !$fasta_seq){
	       $vdj = $toks[$names{"sequence_alignment_aa_heavy"}]."\n".$toks[$names{"sequence_alignment_aa_light"}];
	   }

	   # Note that any field can be specified for this fasta file using seqFasta argument
	   if ($vdj eq ""){
	       $vdj = $toks[$names{$fasta_seq}];
	   }
	   
	   printf FAS ">%s %s\n%s\n", $seq_id, $bfile, $vdj;

       }
   }


}
close(TMP);
close(FS);

sub convertToAA{
    my $ntseq = $_[0];
    my $aaseq = "";

    # Assume in frame
    for (my $c = 0; $c < length($ntseq);$c+=3){
	$aaseq .= $codons_all{substr($ntseq,$c,3)};
    }
    
    return($aaseq);
}
sub readFasta{
    my $fasta = $_[0];
    my %fasta_hash;
    my $header = "";
    my $seq = "";
    open(TMP, "$fasta") or die "Couldn't open $fasta: $?\n";
    while ($line = <TMP>){
	$line =~ s/\n//;
	if ($line =~ /\>(\S+)/){

	    if ($header){
		$fasta_hash{$header} = $seq;
#		printf("%-50s %s %6d %4d\n",$header, $seq);
	    }
	    $seq = "";
	    $header = $1;
	} else {
	    $seq .= $line;
	}
    }
    close(TMP);
    $fasta_hash{$header} = $seq;

    return (%fasta_hash);
}

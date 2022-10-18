#!/usr/bin/env perl

use strict;
use warnings;

#-----------------------------------------------------------------------------
# Global variables
#
my @ARGS = @ARGV;					# Array of passed options
my $VERSION = "";					# Version of this script

#-----------------------------------------------------------------------------
# Load aditional packages
#
use Getopt::Long;					# To easily retrieve arguments from command-line
use Term::ANSIColor qw(:constants);	# Colored output for the terminal

use Data::Dumper;

#-----------------------------------------------------------------------------
# Log and debug
#
sub say 	{ print @_, "\n"; }
sub ohai	{ say BOLD, BLUE,	"==> ", RESET, @_;}
sub oops	{ say BOLD, YELLOW,	"/!\\ ",RESET, @_;}
sub fail	{ say BOLD, RED,	"!!! ", RESET, @_;}

#-----------------------------------------------------------------------------
# Global variables
#
my %OPTIONS;						# Hash of passed options
my $OPTIONS_VALID;					# Are the passed options valid ?
# my @ARGS = @ARGV;					# Array of passed options
my $COMMAND = `basename $0`;		# Name of the script
chomp($COMMAND);
# my $VERSION = "0.1";					# Version of this script


#
# Generic variables
#
my $VERBOSE = 0;					# Verbose mode
my $DEBUG = 0;						# Debug mode

#
# Script specific variables
#
#my $inputFA = $ARGV[0];
#my $inputSAM = $ARGV[1];
#my $maxS = $ARGV[2];

my $sam_fn;
my $fa_fn;
my $maxS;

#-----------------------------------------------------------------------------
# Get passed arguments and check for validity.
#
$OPTIONS_VALID = GetOptions(
	\%OPTIONS,
	'help|h'     => sub { USAGE(); },
	'version'    => sub { VERSION_MESSAGE(); },
	'verbose|v+' => \$VERBOSE,
	'debug|d'    => sub { $VERBOSE = 1; $DEBUG = 1; },
	'sam=s'      => \$sam_fn,
	'fa=s'       => \$fa_fn,
	'maxS=i'     => \$maxS,
);

unless ($OPTIONS_VALID) {
	fail "Error in arguments.";
	USAGE(1);
}

sub USAGE {
	my $exitval = defined($_) ? $_ : 0;
	print <<EOF;
$COMMAND [-h|--help] [--version]
	--help, -h          : Print this help, then exit
	--version           : Print the script version, then exit
EOF
	exit $exitval;
}

sub VERSION_MESSAGE {
	say "This is ", BOLD, "$COMMAND", RESET, " v$VERSION.";
	print <<EOF;
Copyright (c) 2010 Guillaume-Jean Herbiet  (http://herbiet.net)
This is free software; see the source for copying conditions. There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
EOF
	exit 0;
}

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Core code
#
print("start...\n");
{
	my $pick_data = proc_sam_file();
	gen_PS($pick_data);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Additional functions
#
sub gen_PS {
	my $pick_data = shift;
	# ohai Dumper($pick_data);


	# basic value set 
	#my $fontSize = 14;
	my $font_size  = 20;
	my $name_dis_x = 10;
	my $seq_dis_x  = 15;
	my $dis_y	   = 15;
	my $group_dis  = 10;

	my $ps_dir = "psFile";
	`mkdir $ps_dir` if (!-d $ps_dir);
	`rm -f $ps_dir/*`;

	my ($ref_name, $ref_seq) = single_ctg_fasta($fa_fn);
	# ohai $ref_name;
	# ohai $ref_seq;
	my $trans_data = re_position($pick_data);
	# my @transData = @{&TransPositionBySEGA($pick_data)};
	
	my $group_data = group_data($trans_data, $group_dis);
	# my @groupData = @{&GroupData($trans_data, $group_dis)};

	# output file > main table, postscript
	# for my $group_ref (@{$group_data}) {
		# my ($grp_start, $grp_end, $name_maxlen, $sorted_group) = sortByPos($group_ref);
	for my $i (0 .. $#{$group_data}) {
		my ($grp_start, $grp_end, $name_maxlen, $sorted_group) = sort_by_pos($group_data->[$i]);
		
		my $ctg_len = $grp_end - $grp_start;
		my $grp_ctg = join("\t", "$ref_name\_($grp_start..$grp_end)", $ref_name, $grp_start, $ctg_len."M", substr($ref_seq, $grp_start-1, $ctg_len));
		# ohai Dumper $sorted_group;
		unshift @{$sorted_group}, $grp_ctg;
		# ohai Dumper $sorted_group;

		# count start X, Y 
		my $char_x = 50 + ($name_dis_x * $name_maxlen); # name always start from 0 & dis with seq 50
		my $char_y = 0;
		my $max_x  = $char_x + ($grp_end - $grp_start) * $seq_dis_x;
		my $max_y  = $char_y + (scalar @{$sorted_group}) * $dis_y;

		# create ps 
		open(my $outPS, ">", "$ps_dir/output$i.ps")
			or die "Can't open > $ps_dir/output$i.ps: $!";
		my $tagCode = create_position_tag($char_x, $max_y, $seq_dis_x, $dis_y, $grp_start, $grp_end);
		my ($valMatchRef, $valueCode) = create_value_code($char_x, $char_y, $name_dis_x, $seq_dis_x, $dis_y, $grp_end, $sorted_group);		
		my $basicCode = create_basic_code($font_size, $valMatchRef);
		
		print $outPS $basicCode;		
		print $outPS $valueCode;		
		print $outPS $tagCode;
		print $outPS "/psSize [[0 $char_y][$max_x $max_y]] def" ."\n\n";
		print $outPS "init" ."\n";		
		print $outPS "drawChar" ."\n";
		print $outPS "showpage" ."\n";
		close($outPS);
	}
	`cat $ps_dir/*.ps > merge.ps`;
	`ps2pdf merge.ps`; 
}

sub create_basic_code {
	
	my $fontSize =   $_[0];
	my %valMatch = %{$_[1]};

	my $basicCode = '';

	$basicCode .= "/fsize  $fontSize def" ."\n";
	$basicCode .= "/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def" ."\n\n";
	
	$basicCode .= "/init {" ."\n";	
	$basicCode .= "\t". "72 216 translate" ."\n";	
	$basicCode .= "\t". "/xmax -1000 def /xmin 10000 def" ."\n";
	$basicCode .= "\t". "/ymax -1000 def /ymin 10000 def" ."\n\n";
	$basicCode .= "\t". "psSize {" ."\n";
	$basicCode .= "\t\t". "aload pop" ."\n";
	$basicCode .= "\t\t". "dup ymin lt {dup /ymin exch def} if" ."\n";
	$basicCode .= "\t\t". "dup ymax gt {/ymax exch def} {pop} ifelse" ."\n";
	$basicCode .= "\t\t". "dup xmin lt {dup /xmin exch def} if" ."\n";
	$basicCode .= "\t\t". "dup xmax gt {/xmax exch def} {pop} ifelse" ."\n";
	$basicCode .= "\t". "} forall" ."\n\n";	
	$basicCode .= "\t". "/size {xmax xmin sub ymax ymin sub max} bind def" ."\n";
	$basicCode .= "\t". "72 6 mul size div dup scale" ."\n";
	$basicCode .= "\t". "size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div translate" ."\n";	
	$basicCode .= "} bind def" ."\n\n";

	
	# set value function
	my $Mvalue = '';
	my $value  = '';
	foreach my $posName (keys %valMatch){
		my $code = "\t". "[] 0 $posName { aload pop moveto dup $valMatch{$posName} exch 1 getinterval cshow 1 add } forall pop";	
		
		if($posName =~ /M/){ $Mvalue .= "$code\n";
		}else{ $value .= "$code\n";}
	}
	
	
	$basicCode .= "/drawChar {" ."\n";
	$basicCode .= "\t". "/CourierNew findfont fsize scalefont setfont 0 setgray" ."\n";
	$basicCode .= "\t". "[] 0 tagPos { aload pop moveto dup tagStr exch 1 getinterval cshow 1 add } forall pop" ."\n";
	$basicCode .= "\t". "[] 0 intPos { aload pop moveto dup intStr exch 1 getinterval cshow 1 add } forall pop" ."\n";	
	$basicCode .= $value ."\n";
	$basicCode .= "\t". "/CourierNew findfont fsize scalefont setfont 1 0 0 setrgbcolor" ."\n";
	$basicCode .= $Mvalue ."\n";
	$basicCode .= "} bind def" ."\n\n";
	
	return $basicCode;
}

sub create_value_code {
	my ($startX, $startY, $nameDisX, $seqDisX, $disY, $min, $data) = @_;

	my %valSet;
	my %valMatch;
	my $psCode = '';
	my $maxY   = $startY + ($#{$data} * $disY);
		
	for my $i (0 .. $#{$data}) {
		my ($qName, $sName, $start, $sega, $seq) = split(/\t/, $data->[$i]);
		
		my $x = $startX + (($start - $min) * $seqDisX);
		my $y = $maxY - ($i * $disY);
		
		# set name part		
		$valSet{"Name$i"} = $qName;
		$valSet{"NamePos$i"} = create_pos_array(0,  $y, $nameDisX, $qName);
		$valMatch{"NamePos$i"} = "Name$i";
		
		# set position, seqeunce
		my ($valMatchRef, $valSetRef) = classify_position_by_CIGAR($i, $x, $y, $seqDisX, $sega, $seq, \%valMatch, \%valSet);
		%valMatch = %{$valMatchRef};
		%valSet   = %{$valSetRef};
	}
	
	# set value code
	foreach my $name (keys %valSet){
		if($valSet{$name} =~ /^\S+$/){
			$psCode = $psCode ."/$name ($valSet{$name}) def\n";
		}else{
			$psCode = $psCode ."/$name [$valSet{$name}] def\n";
		}
	}
	
	
	return \%valMatch, $psCode;
}

sub classify_position_by_CIGAR {

	my $dataCT	  =   $_[0];
	my $x		  =   $_[1];
	my $y		  =   $_[2];
	my $dis		  =   $_[3];
	my $sega	  =   $_[4];
	my $seq		  =   $_[5];
	my %valMatch  = %{$_[6]};
	my %valSet	  = %{$_[7]};
	
		
	my %seqVal;
	my %posVal;
	my $val   = '';
	my $index = 0;
	my $len   = 0;
	
	for my $i (0..length($sega)-1){
		my $char = substr($sega, $i, 1);
		
		if($char =~ /\d/){
			$val = $val . $char;
		
		}else{					
			$index = $index + $len;
			$len   = ($char eq "D")? 0 : $val;
			$val   = '';
			
			if($char ne "D"){
				my $nowX   = $x + ($index* $dis);
				my $tagSeq = substr($seq, $index, $len);				
				my $array  = create_pos_array($nowX, $y, $dis, $tagSeq);
				
				$seqVal{$char} = (defined($seqVal{$char}))? $seqVal{$char}.$tagSeq : $tagSeq;
				$posVal{$char} = (defined($posVal{$char}))? $posVal{$char}.$array  : $array;
			}
		}
	}
	
	# only pick M case	
	foreach my $key (keys %seqVal){
		if($key eq 'M'){
			$valSet{"SeqM$dataCT"} = (defined($valSet{"SeqM$dataCT"}))? $valSet{"SeqM$dataCT"}.$seqVal{$key} : $seqVal{$key};
			$valSet{"PosM$dataCT"} = (defined($valSet{"PosM$dataCT"}))? $valSet{"PosM$dataCT"}.$posVal{$key} : $posVal{$key};			
			$valMatch{"PosM$dataCT"} = "SeqM$dataCT";
			
		}else{
			$valSet{"Seq$dataCT"} = (defined($valSet{"Seq$dataCT"}))? $valSet{"Seq$dataCT"}.$seqVal{$key} : $seqVal{$key};
			$valSet{"Pos$dataCT"} = (defined($valSet{"Pos$dataCT"}))? $valSet{"Pos$dataCT"}.$posVal{$key} : $posVal{$key};			
			$valMatch{"Pos$dataCT"} = "Seq$dataCT";
		}
	}
		
	return \%valMatch, \%valSet;
}

sub create_pos_array {
	my ($startX, $startY, $dis, $sequence) = @_;

	my $arrayStr = '';
		
	for my $i (0 .. length($sequence)-1){	
		my $nowX = $startX + ($i*$dis);		
		$arrayStr .= "[$nowX $startY]";	
	}
	
	return $arrayStr;
}

sub create_position_tag {
	my ($x, $maxY, $disX, $disY, $min, $max)  = @_;

	my @tagStr;
	my @tagPos;
	my $nowX   = 0;
	my $ct     = 0;
	
	$tagStr[0] = $min;
	$tagPos[0] = $x;
	for my $i($min+1 .. $max-1) {
		$nowX = $x+ (($i-$min) * $disX); $ct++;
		if ($ct == 100) {
			push @tagStr, $i;
			push @tagPos, $nowX;
			$ct=0;
		}
	}
	
	if ($tagPos[$#tagPos] != $nowX) {
		$tagStr[$#tagStr+1] = $max;
		$tagPos[$#tagPos+1] = $nowX;
	}
	
	my $tagCode = '';
	
	$tagCode .= "/tagStr (". '|'x($#tagStr+1) .") def" ."\n";
	$tagCode .= "/tagPos [";
	for my $i (0 .. $#tagPos) {
		$tagCode .= "[". $tagPos[$i] ." ". $maxY ."]";
	}
	$tagCode .= "] def" ."\n";

	my $intStr ='';
	my $intPos ='';
	my $intDis = $disX-5;
	for my $i (0..$#tagPos) {
		$intStr .= $tagStr[$i];
		foreach my $j(0..length($tagStr[$i])-1){
			$intPos .= "[". ($tagPos[$i]+($j*$intDis)) ." ". ($maxY+$disY) ."]";
		}
	}
	
	$tagCode .= "/intStr ($intStr) def" ."\n";
	$tagCode .= "/intPos [$intPos] def" ."\n";
	
	return $tagCode;
}

#
# modified from the original SortBySEGApos()
# 1. sort the contig by start position
# 2. get the first position of a group
# 3. get the last position of a group
# 4. get the longest length of the contig name
#
sub sort_by_pos {
	my $group_ref = shift;

	# ohai Dumper($group_ref);

	my @pick_fields = qw(QNAME RNAME POS CIGAR SEQ);

	my @sorted_group = 
		map {${$_->[1]}}
		sort {$a->[0] <=> $b->[0] || ${$a->[1]} cmp ${$b->[1]}}
		map {[(split /\t/)[2], \$_]} @{$group_ref};
		# ohai Dumper(@sorted_group);

 	my %al;
	@al{@pick_fields} = split /\t/, $sorted_group[0];
	my $group_start = $al{POS};
	my $group_end = $group_start + length($al{SEQ});;
	my $name_maxlen = length($al{QNAME});

	for my $entry (@sorted_group) {
		my %align;
		@align{@pick_fields} = split(/\t/, $entry);
		my $start = $align{POS};
		my $end = $start + length($align{SEQ});
		my $name_len = length($align{QNAME});
		$group_end = $end if ($end > $group_end);
		$name_maxlen = $name_len if ($name_len > $name_maxlen);
	}
	return ($group_start, $group_end, $name_maxlen, \@sorted_group);
}


#
# grouping the reads
#
sub group_data {
	my $trans_data = shift;
	my $group_dis = shift;

	my @pick_fields = qw(QNAME RNAME POS CIGAR SEQ);

	my %hit_map;
	for my $entry (@{$trans_data}) {
		my %align;
		@align{@pick_fields} = split(/\t/, $entry);
		my $start = $align{POS};
		my $end = $start + length($align{SEQ}) - 1;

		# ohai "$start, $end";
		
		for my $i ($start .. $end) {
			$hit_map{$i}  = 1;
		}
	}
	
	# my @all_pos = sort {$a <=> $b} keys %hit_map;
	# ohai join " ", @all_pos;
	
	my @group;
	my $pre_pos;
	for my $pos (sort {$a <=> $b} keys %hit_map) {
		# ohai $key;
		if (defined($pre_pos) && $pos - $pre_pos > $group_dis) {
		 	push(@group, $pre_pos);
		}
		$pre_pos = $pos;
	}
	if ($group[$#group] != $pre_pos) {
		push(@group, $pre_pos);
	}
	# ohai Dumper @group;

	my @group_data;
	for my $entry (@{$trans_data}) {
		my %align;
		@align{@pick_fields} = split(/\t/, $entry);
		my $start = $align{POS};
		my $end = $start + length($align{SEQ}) - 1;

		for my $j (0 .. $#group) {
			if ($end <= $group[$j]) {
				# ohai "$end <= $group[$j] $entry";
				push(@{$group_data[$j]}, $entry);
				last;
			}
		}
	}
	# say Dumper @group_data;
	
	return \@group_data;
}

sub re_position {
	my $pick_data = shift;
	# ohai Dumper $pick_data;
	
	my @pick_fields = qw(QNAME RNAME POS CIGAR SEQ);

	my @data;
	for my $entry (@{$pick_data}) {
		# ohai $entry;
		my %align;
		@align{@pick_fields} = split(/\t/, $entry);
		# ohai $align{CIGAR};

		# CIGAR regex \*|([0-9]+[MIDNSHPX=])+
		# page 6 in https://samtools.github.io/hts-specs/SAMv1.pdf
		# page 8 CIGAR operations: ops that consumes {query, reference} = {yes, no}
		my (@ops) = ($align{CIGAR} =~ m/(\*|[0-9]+[MIDNSHPX=])/g);
		my $pos = $align{POS};
		my $start = $pos;
		{
			# only check the first CIGAR operation
			$ops[0] =~ /([0-9]+)(.)/;
			my ($len, $op) = ($1, $2);
			# ohai "$1, $2";
			if ($op =~ /[SI]/) {
				$start = $pos - $len;
			}
		}
		# ohai "$pos -> $start : ", ($pos - $start), " ", $align{CIGAR};

		$align{POS} = $start;
		push @data, join "\t", @align{@pick_fields};
	}
	return \@data;
}

sub single_ctg_fasta {
	my $fn = shift;

	my $ctg_name;
	my $ctg_seq;

	open(my $fh, "<", $fn)
		or die "Can't open < $fn: $!";
	{
		local $/ = undef;
		my $text = <$fh>;
		# ohai $text;
		my @ctgs = split(/>/, $text);
		my @tmp = split(/\s/, $ctgs[1]);
		$ctg_name = shift @tmp;
		$ctg_seq = join '', @tmp;
		# ohai $ctg_name;
		# ohai $ctg_seq;
	}
	close($fh);
	return $ctg_name, $ctg_seq;
}

sub proc_sam_file {
	my $count = 0;
	my @pick_data;

	# output max S part sequence ( S must > $maxS )
	open(my $trimFa, ">", "trim_S.fa") 
		or die "Can't open > trim_S.fa: $!";
	# output pick info table
	open(my $trimInfo, ">", "trim_info.txt") 
		or die "Can't open > trim_info.fa: $!";
	# print info header
	print $trimInfo "Query\tHit\tFLAG\tPosition\tCIGAR\tQuerySequence\n";
	
	say($sam_fn, "\n");
	open(my $samFile, "<", $sam_fn)
		or die "Can't open < $sam_fn: $!";

	# ref: https://en.wikipedia.org/wiki/SAM_(file_format)
	# ref: https://samtools.github.io/hts-specs/SAMv1.pdf
	# ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/
	my @sam_fields      = qw(QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL);
	my @trimInfo_fields = qw(QNAME RNAME FLAG POS CIGAR SEQ);
	my @pick_fields     = qw(QNAME RNAME POS CIGAR SEQ);

	while (<$samFile>) {
		chomp;
		# say ++$count, "\n", $_,;
		# trim header first
		next if ($_ =~ /^@/);
				
		my %align;
		@align{@sam_fields} = split(/\t/, $_);
		
		# ohai join(" ", @align{@sam_fields});
		my ($Sindex, $Slen) = proc_cigar($align{CIGAR});
		# ohai $Sindex, " ", $Slen;

		if($Slen > $maxS) {
			# my $tmp = join " ", @align{@trimInfo_fields};
			# oops $tmp;
			# oops join "\t", @align{@trimInfo_fields}, "\n";
			# oops ">$align{QNAME}\_$Slen", "S\n", substr($align{SEQ}, $Sindex, $Slen), "\n";	
			print $trimInfo join "\t", @align{@trimInfo_fields}, "\n";
			print $trimFa ">$align{QNAME}\_$Slen","S\n". substr($align{SEQ}, $Sindex, $Slen) ."\n";	
			
			# push (@pick_data, "$qName\t$sName\t$pos\t$SEGA\t$qSeq");
			# ohai join "\t", @align{@pick_fields};
			push @pick_data, join "\t", @align{@pick_fields};
		}
	}
	close($samFile);
	close($trimFa);
	close($trimInfo);

	return \@pick_data;
}
sub proc_cigar {
	my $cigar = shift;
	# ohai $cigar;

	# CIGAR regex \*|([0-9]+[MIDNSHPX=])+
	# page 6 in https://samtools.github.io/hts-specs/SAMv1.pdf
	# page 8 CIGAR operations: ops that consumes queries
	my (@ops) = ($cigar =~ m/(\*|[0-9]+[MIDNSHPX=])/g);
	# ohai join " ", @ops;
	my (@query_ops) = grep {/[MISX=]$/} @ops;
	# ohai join " ", @query_ops;

	my $pMaxS = 0;
	my $maxS = 0;
	my $pos = 0;
	for my $op (@query_ops) {
		$op =~ /([0-9]+)(.)/;
		if ($2 eq "S") {
			($pMaxS, $maxS) = ($pos, $1) if ($1 > $maxS);
		}
		$pos += $1;
	}
	# ohai "$pMaxS $maxS";
	return ($pMaxS, $maxS);
}
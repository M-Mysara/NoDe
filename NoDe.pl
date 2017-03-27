#! /usr/local/bin/perl

use strict;
open FH,$ARGV[0];
my @fasta=<FH>;
close FH;
open FH,$ARGV[1];
my @flow=<FH>;
close FH;
open FH,$ARGV[2];
my @flow_ind=<FH>;
close FH;
open FH,$ARGV[3];
my @qual=<FH>;
close FH;
open FH,$ARGV[4];
my @names=<FH>;
close FH;
open FH1,$ARGV[5];
my @align=<FH1>;
close FH1;
open FH1,$ARGV[6];
my @name=<FH1>;
close FH1;
my $oo=$ARGV[7];
my $ID="";
my $i=13;
open FH, ">./temp_split/$oo.Test.arff";
print FH '@relation sirnahabal

@attribute Pre_pre_flow numeric
@attribute Pre_flow numeric
@attribute flow numeric
@attribute Pos_flow numeric
@attribute Pos_pos_flow numeric
@attribute Position numeric
@attribute Phred_before numeric
@attribute Phred numeric
@attribute Phred_after numeric
@attribute PreHomo {N,A,B,C,D,E,F,G,Z}
@attribute Homo {N,A,B,C,D,E,F,G,Z}
@attribute PostHomo {N,A,B,C,D,E,F,G,Z}
@attribute CFE {F,T}
@attribute Final_Results {T,N}

@data
';
close FH;
open FH, ">>./temp_split/$oo.Test.arff";
my @seq_start;
my @seq_dir;
my @orig_ID;
my @orig_Seq;
my @Main_counter;
my $starter=0;
for(my $o=0;$o<scalar(@name);$o++){
	my $freq=scalar(split(",",$name[$o]));
	if(0){}
	else{
		$ID=$name[$o];
		#print $ID."\n";
		$ID=~/([\w]+)	/;
		$ID=$1;
		my $qflag=0;
		for($i=0;$i<scalar(@names)&&$qflag==0;$i=$i+1){
			if($names[$i]=~/$ID/){
				my $temp_fasta=$fasta[$i];
			    push(@orig_ID,$align[($o*2)]);
	 			push(@orig_Seq,$align[($o*2)+1]);
				push(@Main_counter,0);
				#####
				chomp($fasta[$i]);
				my $Seq=$fasta[$i];
				$Seq=substr($Seq,7);
				my $start=0;
				my $end=0;
				my $seq_=$align[($o*2)+1];
				chomp($seq_);
				$seq_=~s/\.//g;
				$seq_=~s/\-//g;
				$starter=length($seq_);
				push(@seq_start,$starter);
				my $marker=0;
				my $start_Seq="";
				my $seq_temp=$seq_;				
				if ($seq_temp=~/^([A,G,T,C]+)n/i){$seq_temp=$1;}
				if($Seq=~/^([\w]+)$seq_temp/i){$start=length($1);$start_Seq=$1;}
				else{
					$marker=1;
					$seq_temp=reverse($seq_temp);
					$seq_temp=~s/T/F/gi;
					$seq_temp=~s/A/T/gi;
					$seq_temp=~s/F/A/gi;
					$seq_temp=~s/C/F/gi;
					$seq_temp=~s/G/C/gi;
					$seq_temp=~s/F/G/gi;
					if($Seq=~/([\w]+)$seq_temp/i){$start=length($1);$start_Seq=$1;}
					else{print "Something is wrong in the sff $ID $seq_	$Seq $i $o \n";exit;}
				}
				push(@seq_dir,$marker);
				$end=$start+length($seq_);
				my @Seq_Arr=split("",$Seq);
				#####
				chomp($flow[$i]);
				my $Flow=$flow[$i];
				$Flow=substr($Flow,10);	
				my @Flow_Arr=split("	",$Flow);
				my $Flow_Index=$flow_ind[$i];
				$Flow_Index=substr($Flow_Index,14);
				my @Flow_Index_Arr=split("	",$Flow_Index);
			#	#####
				chomp($qual[$i]);
				my $Phred_Sc=$qual[$i];
				$Phred_Sc=substr($Phred_Sc,16);
				my @Phred_Sc_Arr=split("	",$Phred_Sc);
			#	
				########
				my $Homo=0;#my $buffer=0;#
				my $CFE_A="";my $CFE_T="";my $CFE_G="";my $CFE_C=""; #carry forward events
				my $CFE_results="T";
				my $Pre_Homo="";
				my $Pos_Homo="";
				my $new_length=0;
				my $Homo_length=0;
				my $Homo_prev_length=0;
				for(my $j=$start;$j<$end;$j++){
					if($j==$start){$Pre_Homo=0;if($Seq_Arr[$j]eq$Seq_Arr[$j+1]){$Homo=1;}else{$Homo=0;}if($Seq_Arr[$j+1]eq$Seq_Arr[$j+2]){$Pos_Homo=$Homo+1;}else{$Pos_Homo=$Homo;}}
					else{
						$Pre_Homo=$Homo;$Homo=$Pos_Homo;
						if($j==$end-1){$Pos_Homo=0;}
						elsif($Seq_Arr[$j+1]eq$Seq_Arr[$j+2]){if($Seq_Arr[$j]eq$Seq_Arr[$j+1]){$Pos_Homo=$Homo+1;}else{$Pos_Homo=1;}}#3shan lessa e7na 3mleen $Pre_Homo=$Homo;
						else{if($Seq_Arr[$j]eq$Seq_Arr[$j+1]){$Pos_Homo=$Homo+1;}else{$Pos_Homo=0;}}
					}
					

					if($Homo>0){$CFE_results="F";if($Seq_Arr[$j]=~/A/i){$CFE_A="T";}if($Seq_Arr[$j]=~/T/i){$CFE_T="T";}if($Seq_Arr[$j]=~/G/i){$CFE_G="T";}if($Seq_Arr[$j]=~/C/i){$CFE_C="T";}}
					else{
						if($Seq_Arr[$j]=~/A/i){if($CFE_A eq "T"){$CFE_results="T";$CFE_A="F";}else{$CFE_results="F";$CFE_A="F";}}
						elsif($Seq_Arr[$j]=~/G/i){if($CFE_G eq "T"){$CFE_results="T";$CFE_G="F";}else{$CFE_results="F";$CFE_G="F";}}
						elsif($Seq_Arr[$j]=~/C/i){if($CFE_C eq "T"){$CFE_results="T";$CFE_C="F";}else{$CFE_results="F";$CFE_C="F";}}
						elsif($Seq_Arr[$j]=~/T/i){if($CFE_T eq "T"){$CFE_results="T";$CFE_T="F";}else{$CFE_results="F";$CFE_T="F";}}
					}
#				}
					my $Count=$j;
					if($Phred_Sc_Arr[$j+1]){}else{$Phred_Sc_Arr[$j+1]=0;}
				my $curr_flow; my $pre_flow=0;my $pre_pre_flow=0;my $post_flow=0;my $post_post_flow=0;
					my $P_1=0;my $P_2=0;my $P_3=0;my $N_1=0;my $N_2=0;my $N_3=0;	
					if(($Flow_Index_Arr[$j-1]-1)==($Flow_Index_Arr[$j]-2)){$pre_flow=$Flow_Arr[$Flow_Index_Arr[$j]-2];$pre_pre_flow=0;
					$N_3=0;$N_2=0;$N_1=$Flow_Arr[$Flow_Index_Arr[$j]-2];}
				 elsif(($Flow_Index_Arr[$j-1]-1)==($Flow_Index_Arr[$j]-3)){$pre_flow=$Flow_Arr[$Flow_Index_Arr[$j]-3];$pre_pre_flow=$Flow_Arr[$Flow_Index_Arr[$j]-2];
				$N_3=0;$N_2=$Flow_Arr[$Flow_Index_Arr[$j]-3];$N_1=$Flow_Arr[$Flow_Index_Arr[$j]-2];}
				 elsif(($Flow_Index_Arr[$j-1]-1)==($Flow_Index_Arr[$j]-4)){$pre_flow=$Flow_Arr[$Flow_Index_Arr[$j]-4];if($Flow_Arr[$Flow_Index_Arr[$j]-2]>$Flow_Arr[$Flow_Index_Arr[$j]-3]){$pre_pre_flow=$Flow_Arr[$Flow_Index_Arr[$j]-2];}else{$pre_pre_flow=$Flow_Arr[$Flow_Index_Arr[$j]-3];}
$N_3=$Flow_Arr[$Flow_Index_Arr[$j]-4];$N_2=$Flow_Arr[$Flow_Index_Arr[$j]-3];$N_1=$Flow_Arr[$Flow_Index_Arr[$j]-2];}
				    if(($Flow_Index_Arr[$j+1]-1)==($Flow_Index_Arr[$j])){$post_flow=$Flow_Arr[$Flow_Index_Arr[$j]];$post_post_flow=0;
				$P_3=0;$P_2=0;$P_1=$Flow_Arr[$Flow_Index_Arr[$j]];}
				 elsif(($Flow_Index_Arr[$j+1]-1)==($Flow_Index_Arr[$j]+1)){$post_flow=$Flow_Arr[$Flow_Index_Arr[$j]+1];$post_post_flow=$Flow_Arr[$Flow_Index_Arr[$j]];
				$P_3=0;$P_2=$Flow_Arr[$Flow_Index_Arr[$j]+1];$P_1=$Flow_Arr[$Flow_Index_Arr[$j]];}	
			     elsif(($Flow_Index_Arr[$j+1]-1)==($Flow_Index_Arr[$j]+2)){$post_flow=$Flow_Arr[$Flow_Index_Arr[$j]+2];if($Flow_Arr[$Flow_Index_Arr[$j]+1]>$Flow_Arr[$Flow_Index_Arr[$j]]){$post_post_flow=$Flow_Arr[$Flow_Index_Arr[$j]+1];}else{$post_post_flow=$Flow_Arr[$Flow_Index_Arr[$j]];}
$P_3=$Flow_Arr[$Flow_Index_Arr[$j]+2];$P_2=$Flow_Arr[$Flow_Index_Arr[$j]+1];$P_1=$Flow_Arr[$Flow_Index_Arr[$j]];}
				my $Homo_="";my $Pos_Homo_="";my $Pre_Homo_="";
				if($Homo eq 0){$Homo_="N";}elsif($Homo eq 1){$Homo_="A";}elsif($Homo eq 2){$Homo_="B";}elsif($Homo eq 3){$Homo_="C";}elsif($Homo eq 4){$Homo_="D";}elsif($Homo eq 5){$Homo_="E";}elsif($Homo eq 6){$Homo_="F";}elsif($Homo eq 7){$Homo_="G";}
				if($Pos_Homo eq 0){$Pos_Homo_="N";}elsif($Pos_Homo eq 1){$Pos_Homo_="A";}elsif($Pos_Homo eq 2){$Pos_Homo_="B";}elsif($Pos_Homo eq 3){$Pos_Homo_="C";}elsif($Pos_Homo eq 4){$Pos_Homo_="D";}elsif($Pos_Homo eq 5){$Pos_Homo_="E";}elsif($Pos_Homo eq 6){$Pos_Homo_="F";}elsif($Pos_Homo eq 7){$Pos_Homo_="G";}
				if($Pre_Homo eq 0){$Pre_Homo_="N";}elsif($Pre_Homo eq 1){$Pre_Homo_="A";}elsif($Pre_Homo eq 2){$Pre_Homo_="B";}elsif($Pre_Homo eq 3){$Pre_Homo_="C";}elsif($Pre_Homo eq 4){$Pre_Homo_="D";}elsif($Pre_Homo eq 5){$Pre_Homo_="E";}elsif($Pre_Homo eq 6){$Pre_Homo_="F";}elsif($Pre_Homo eq 7){$Pre_Homo_="G";}
				if($Homo>$Pos_Homo){if($Homo_ eq "N"){}else{$Homo_="Z";}}
				if($Pre_Homo>$Homo){if($Pre_Homo_ eq "N"){}else{$Pre_Homo_="Z";}}
				if($Pos_Homo=~/[0,1]+/){}else{if($Seq_Arr[$j+1]eq $Seq_Arr[$j+2]){}else{$Pos_Homo_="Z";}}
				my $prev_position=$Seq_Arr[$j-1];
				if($Homo != 0){
					if($Seq_Arr[$j]eq $Seq_Arr[$j+1]){
						if($j eq $start){if($Seq_Arr[$j]eq $prev_position){$Homo_length++;$prev_position="N";}}
						if($Seq_Arr[$j]eq $prev_position){$Homo_length=$Homo_length;}
						else{
							$Homo_length=$Homo_length+2;
							if($Seq_Arr[$j+1]eq $Seq_Arr[$j+2]){
								$Homo_length++;
								if($Seq_Arr[$j+2]eq $Seq_Arr[$j+3]){
									$Homo_length++;
									if($Seq_Arr[$j+3]eq $Seq_Arr[$j+4]){
										$Homo_length++;
										if($Seq_Arr[$j+4]eq $Seq_Arr[$j+5]){
											$Homo_length++;
											if($Seq_Arr[$j+5]eq $Seq_Arr[$j+6]){
												$Homo_length++;
											}else{$Homo_length=$Homo_length;}
										}else{$Homo_length=$Homo_length;}
									}else{$Homo_length=$Homo_length;}
								}else{$Homo_length=$Homo_length;}
							}else{$Homo_length=$Homo_length;}
						}
					}
					elsif($Seq_Arr[$j]eq $Seq_Arr[$j-1]){$Homo_length=$Homo_length;}
				}
				else{
					$Homo_length=0;
				}
                                my $Post_Homo_length=0;
				if($Pos_Homo == 1){
					$Post_Homo_length=$Post_Homo_length+2;
					if($Seq_Arr[$j+2]eq $Seq_Arr[$j+3]){
						$Post_Homo_length++;
						if($Seq_Arr[$j+3]eq $Seq_Arr[$j+4]){
							$Post_Homo_length++;
							if($Seq_Arr[$j+4]eq $Seq_Arr[$j+5]){
								$Post_Homo_length++;
								if($Seq_Arr[$j+5]eq $Seq_Arr[$j+6]){
										$Post_Homo_length++;
									if($Seq_Arr[$j+6]eq $Seq_Arr[$j+7]){
										$Post_Homo_length++;
									}else{$Post_Homo_length=$Post_Homo_length;}
								}else{$Post_Homo_length=$Post_Homo_length;}
							}else{$Post_Homo_length=$Post_Homo_length;}
						}else{$Post_Homo_length=$Post_Homo_length;}												
					}else{$Post_Homo_length=$Post_Homo_length;}
				}
				elsif($Pos_Homo>1){$Post_Homo_length=$Homo_length;}
				else{$Post_Homo_length=0;}
				if($Homo_length==0){$curr_flow=$Flow_Arr[$Flow_Index_Arr[$j]-1];}
				else{$curr_flow=$Flow_Arr[$Flow_Index_Arr[$j]-1]/$Homo_length;}
				if($Homo_prev_length==0){$pre_flow=$pre_flow;}
				else{$pre_flow=$pre_flow / $Homo_prev_length;}
				if($Post_Homo_length==0){$post_flow=$post_flow;}
				else{$post_flow=$post_flow / $Post_Homo_length;}
				if($Homo_ eq "N"){}
				elsif($Homo_ eq "A"){$post_flow=$curr_flow;}
				elsif($Homo_ eq "Z"){$pre_flow=$curr_flow;}
				else{$pre_flow=$curr_flow;$post_flow=$curr_flow;}
				print FH $pre_pre_flow.",",$pre_flow.",".$curr_flow.",".$post_flow.",".$post_post_flow.",";
				$Homo_prev_length=$Homo_length;
				print FH $Count.",".$Phred_Sc_Arr[$j-1].",".$Phred_Sc_Arr[$j].",".$Phred_Sc_Arr[$j+1].",".$Pre_Homo_.",".$Homo_.",".$Pos_Homo_.",".$CFE_results.","."T\n";	
				}
				
		 		$qflag++;
			}
			else{next;}
		}
	}
}
close FH;
my $model='NoDe.model';my $weka='SMO';
#print "java -Xmx2000M -classpath ./weka-3-6-8/weka.jar weka.classifiers.functions.$weka -l ./$model -T ./temp_split/$oo.Test.arff -p 0 \> ./temp_split/$oo.Test.Final";
system"java -Xmx2000M -classpath ./weka.jar weka.classifiers.functions.$weka -l ./$model -T ./temp_split/$oo.Test.arff -p 0 > ./temp_split/$oo.Test.Final";
open WEKA,"./temp_split/$oo.Test.Final";
my $results="";
while(my $weka=<WEKA>){
	if($weka=~/[\s]+\d\:\w[\s]+\d\:(\w)[\s]+/){$results=$results.$1.',';}
}
close WEKA;
chop($results);
my @res_arr=split(",",$results);
for(my $i=0;$i<scalar(@orig_ID);$i++){
	print $orig_ID[$i];
	if($Main_counter[$i]==1){print$orig_Seq[$i];}
	else{
		my $seq=$orig_Seq[$i];
		chomp($seq);
		my @arr_seq=split("",$seq);
		my @temp_seq=splice(@res_arr,0,$seq_start[$i]);
		my $h=0;
		if($seq_dir[$i]==1){@temp_seq=reverse(@temp_seq);}
		for(my $l=0;$l<scalar(@arr_seq);$l++){
			if($arr_seq[$l]eq'.'){print '.';}
			elsif($arr_seq[$l]eq'-'){print '-';}
			elsif($temp_seq[$h]eq'T'){$h++;print$arr_seq[$l];}
			elsif($temp_seq[$h]eq'S'){$h++;print"N";}
			elsif($temp_seq[$h]eq'I'){$h++;print"N";}
			elsif($temp_seq[$h]eq'D'){$h++;print"N";}
			elsif($temp_seq[$h]eq'N'){$h++;print"N";}
		}
		print "\n";
	}
}

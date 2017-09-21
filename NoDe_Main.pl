#! /usr/local/bin/perl
#.....................License.................................
#	 NoDe a Software for NextGeneration sequencing data denoising
#    Copyright (C) 2014  <M.Mysara et al>

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
	
#......................Packages Used..........................
use strict;
use Getopt::Std;
my %opts;
my $diffs;
my $group;
getopt('npfsgd',\%opts);

if(!defined($opts{f}) && !defined($opts{n}) && !defined($opts{s})){usage();exit;}
if(!defined($opts{p})){$opts{p}=1;}
if(!defined($opts{d})){
	open FASTA,$opts{f};
	my @FASTA=<FASTA>;
	close FASATA;
	my $sequence="";
	my $count=0;
	for(my $k=0;$k<scalar(@FASTA);$k++){
		if($FASTA[$k]=~/>/){$count++;}
		else{
			chomp($FASTA[$k]);
			$sequence=$sequence.$FASTA[$k];}
	}
	$sequence=~s/\.//gi;
	$sequence=~s/\-//gi;
	my $length=length($sequence)/$count;
	my $qflage=0;
	for(my $p=100;$qflage==0;$p=$p+100){
		if($length>$p){$diffs=$p/100;}
		else{
			#$diffs=$diffs+1;
			$qflage=1;
		}
	}
}
else{$diffs=$opts{d};}
if(defined($opts{g})){$group=$opts{g};}
open FH,$opts{f};
my @seq=<FH>;
close FH;
open FH,$opts{n};
my @name=<FH>;
close FH;
my $core=$opts{p};
my $logfile=int(rand(1000000000000000));
$logfile='NoDe_'.$logfile.'.logfile';
my $seq_len=scalar(@seq)/2;
my $seq_frac=int($seq_len/$core);
my $reminder=$seq_len%$core;
system "rm -rf temp_split";
mkdir("temp_split");
my $i=0;
my $file='./temp_split/'.$i.'.fasta';
my $file1='./temp_split/'.$i.'.name';
foreach $i(0..$core-1){
	$file='./temp_split/'.$i.'.fasta';
        $file1='./temp_split/'.$i.'.name';
        open FH,">>",$file;
        open FH1,">>",$file1;
        for(my $j=($i*$seq_frac*2);$j<($i+1)*$seq_frac*2;$j++){
        	print FH $seq[$j];
		my $remind=$j%2;
		if($remind){}
		else{print FH1 $name[($j/2)];}
        }
        close FH;
	close FH1;
}
open FH, ">>",$file;
open FH1,">>",$file1;
for(my $j=($core)*$seq_frac*2;$j<($core*$seq_frac*2)+($reminder*2);$j++){
	print FH $seq[$j];
	my $remind=$j%2;
        if($remind){}
        else{print FH1 $name[($j/2)];}

}
close FH;
close FH1;
#####greping from sff file
system " grep -P \"^Bases: \" $opts{s} > ./temp_split/fasta";
system " grep \"Flowgram: \" $opts{s} > ./temp_split/flow";
system " grep \"Flow Indexes: \" $opts{s} > ./temp_split/flow_ind";
system " grep \"Quality Scores: \" $opts{s} > ./temp_split/qual";
system " grep \">\" $opts{s} > ./temp_split/name";
my $fasta='./temp_split/fasta';
my $flow='./temp_split/flow';
my $flow_ind='./temp_split/flow_ind';
my $qual='./temp_split/qual';
my $name='./temp_split/name';
###########
foreach my $child (0..$core-1) {
	my $file_='./temp_split/'.$child.'.Result';
        unlink $file_;# or warn "Could not unlink $file_: $!";
        my $file='./temp_split/'.$child.'.fasta';
	my $file1='./temp_split/'.$child.'.name';
	my $pid = fork();
        if ($pid == -1) {
        	die;
        }
	elsif ($pid == 0) {
		print  "perl NoDe.pl $fasta $flow $flow_ind $qual $name $file $file1 $child \> $file_ \n";
                exec "perl NoDe.pl $fasta $flow $flow_ind $qual $name $file $file1 $child > $file_" or die;                 
        }
}
while (wait() != -1) {};
print "Done\n";
my $cat="";
foreach my $child (0..$core-1) {
	my $file_='./temp_split/'.$child.'.Result';
	$cat=$cat." ".$file_;
}
system "cat $cat > ./temp_split/Results.fasta";
system "cp $opts{n} ./temp_split/Results.names";
#running SLP modified algorithm
if ($opts{g}){
	system "\.\/mothur \"\#pre\.cluster(fasta=./temp_split/Results.fasta,name=./temp_split/Results.names,diffs=$diffs,group=$group,processors=$core\)\">> $logfile";
}
else{
	system "\.\/mothur \"\#pre\.cluster(fasta=./temp_split/Results.fasta,name=./temp_split/Results.names,diffs=$diffs,processors=$core\)\">> $logfile";
}
#Getting reads without the Markers
system "cut -f1 ./temp_split/Results.precluster.names > ./temp_split/Results.precluster.accnos";
system "cp $opts{f} ./temp_split/Results.precluster_.fasta";
system "\.\/mothur \"\#get\.seqs(fasta=./temp_split/Results.precluster_.fasta,accnos=./temp_split/Results.precluster.accnos\)\" >> $logfile";
unlink("./temp_split/Results.precluster_.fasta");
mkdir('NoDe_Final');

system "cp ./temp_split/Results.precluster_.pick.fasta ./NoDe_Final/Results.NoDe.fasta";
system "cp ./temp_split/Results.precluster.names ./NoDe_Final/Results.NoDe.names";

sub usage{
print
"	
	||||||||||||||||||||||||||||||||||||||||||||||||
	||              Welcome To NoDe		      ||
	|| A Software For Sequencing Data Denoising   ||
	||    Copyright (C) 2013  <M.Mysara et al>    ||
	||||||||||||||||||||||||||||||||||||||||||||||||

 NoDe version 1, Copyright (C) 2014, M.Mysara et al
 NoDe comes with ABSOLUTELY NO WARRANTY.
 This is free software, and you are welcome to redistribute it under
 certain conditions; type `notepad copying.txt\' for details.;

 The software also includes \"WEKA-3-6\" and \"mothur\" Both under GNU Copyright
 
 
Command Syntax:
(perl) NoDe_Main.pl {options} 	
 Use the Following Mandatory Options:\n
 -s The sff file
 -n Name file with the redundancy
 -f Aligned sequences file (accept only AGTC bases).
	
 Use the Following Non-Mandatory Options:\n
 -p number of processors, defaul 1
 -g group file \(in case of having different sample-groups\)
 -d differences tolerated (Defualt: will be automatically calculated
	to be below 98% distance)
 
 For Queries about the installaion, type \'less README.txt\'
 For Queries about the Copy rights, type \'less COPYING.txt\'
	
";
}

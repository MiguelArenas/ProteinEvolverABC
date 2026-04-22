#!/usr/bin/env perl 
$newdir=1;  # 1 if create one dir for each protein, 0 if all in the same
$cluster=0; # 1 if run in cluster, 0 in your machine
$List="????"; # file with pdb codes
open(my $fi, '<:encoding(UTF-8)', $List)
    or die "Could not open file '$List' $!";


$Input="Input_Prot_evol.in";


while (my $row = <$fi>) {
    chomp $row;
    @res = split(/\s+/, $row);
    $P=$res[0];
    print "Prot= $P \n";

    $tmp=sprintf("%s.tmp", $P);
    `sed "s/PPPP/"$P"/g" $Input > $tmp`;
# If you want to modify other parameters:
#`sed -i "s/LLL/"$L"/g" $tmp`;

# If you want to make a directory:
    if($newdir){
	$dir=$P;
	`mkdir $dir`;
	cp $tmp $dir;
    }

    if($cluster==0){
	`Prot_Evol $tmp > $P.log`; # run in your machine
    }else{
	# Run in the cluster; you have to log in into medusa31
	$script=sprintf("%s.script", $P); 
	print "Script= ", $script,"\n";
	open(my $fh, '>', $script);
	print $fh "Prot_evol $tmp > $P.log\n";
	close $fh;
	`chmod u+x $script`;
	`qsubmit.pl -q x86_64 -s $script`;
    }
    if($newdir){`cd ..`;}
}

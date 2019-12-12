##get initial vcf from snv-process

open(DATA,"<D:/file/willow/snp/new_criterion/new_group2/groups.txt") or die "can't open data";

my %group;
my %number;
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$group{$temp[0]}=$temp[1];
	$number{$temp[1]}++;
}


@file=glob("D:/file/willow/snp/new_criterion/new_group2/*.csv") or die "can't open glob";

for $filename(@file){
	$filename=~/D:\/file\/willow\/snp\/new_criterion\/new_group2\/(.*).csv/;
	open(DATA,"<$filename") or die "can't open data";
	open(INPUT,">$1.stats") or die "can't open input";
	LABEL:while(<DATA>){
		chomp;
		my %hash=();
		@temp=split(/\s+/);
		next if($temp[-1]=~/FPD/);
		next if($temp[-1]=~/MISSING/);
		next unless(/Confidence/);
	  @sample=split(/;/,$temp[2]);
	  
	  if(/FPS=\((.*)\);/){
	  	@fps=split(/,/,$1);
	  	for $fps(@fps){
	  		push(@sample,$fps);
	  	}
	  }
	  
	  if(/GRPS=(.*?);/){
	    @grps=split(/,/,$1);
	    for $grps(@grps){
	    	$grps=~s/\(//;
	    	$grps=~s/\)//;
	    	push(@sample,$grps);
	    }	
	  }
	  
	  next if($#sample==0);
	  for $sample(@sample){
	  	$hash{$group{$sample}}++;
	  }
	  
	  $branch2=0;$branch5=0;$branch6=0;
	  
	  for $key(sort keys %hash){
	  	if($key=~/branch2/){
	  		$diff=$number{$key}-$hash{$key};
	  		$branch2+=$diff;
	  		next;
	  	}
	  	
	  	if($key=~/branch5/){
	  		$diff=$number{$key}-$hash{$key};
	  		$branch5+=$diff;
	  		next;
	  	}
	  	
	  	if($key=~/branch6/){
	  		$diff=$number{$key}-$hash{$key};
	  		$branch6+=$diff;
	  		next;
	  	}
	  	next LABEL unless($hash{$key} eq $number{$key}); 
	  }
	  
	  next LABEL if($branch2>5);next LABEL if($branch5>8);next LABEL if($branch6>13);##tolerate 50% samples in every branch discard mutation 
	  
	  for $key(sort keys %hash){
	  	print INPUT "$key:$hash{$key}/$number{$key};";
	  }
	  print INPUT "\t$_\n";
	}
	close(DATA);
	close(INPUT);
}


@file=glob("D:/file/willow/snp/new_criterion/new_group2/*.stats") or die "can't open glob";

open(INPUT,">D:/file/willow/snp/new_criterion/new_group2/all.txt");
for $filename(@file){
	open(DATA,"<$filename");
	while(<DATA>){
		chomp;
		print INPUT "$_\n";
	}
	close(DATA);
}
##get initial vcf from snv-process

open(DATA,"<D:/file/willow/newref/snp_newgroup/groups.txt") or die "can't open data";

my %group;
my %number;
while(<DATA>){
	chomp;
	@temp=split(/\s+/);
	$group{$temp[0]}=$temp[1];
	$number{$temp[1]}++;
}


@file=glob("D:/file/willow/newref/snp_newgroup/*.csv") or die "can't open glob";

for $filename(@file){
	$filename=~/D:\/file\/willow\/newref\/snp_newgroup\/(.*).csv/;
	open(DATA,"<$filename") or die "can't open data";
	open(INPUT,">$1.stats") or die "can't open input";
	LABEL:while(<DATA>){
		chomp;
		my %hash=();
		@temp=split(/\s+/);
		next if($temp[-1]=~/FPD/);
		next if($temp[-1]=~/MISSING/);
		next unless(/Confidence/); ###只保留confidence的结果
	  @sample=split(/;/,$temp[2]);
	  
	  if(/FPS=(.*?);/ && !/FPS=NA/){
	  	@fps=split(/,/,$1);
	  	for $fps(@fps){
	  		$fps=~s/\(//;
	    	$fps=~s/\)//;
	  		push(@sample,$fps);
	  	}
	  }
	  
	  if(/GRPS=(.*?);/ && !/GRPS=NA/){
	    @grps=split(/,/,$1);
	    for $grps(@grps){
	    	$grps=~s/\(//;
	    	$grps=~s/\)//;
	    	push(@sample,$grps);
	    }	
	  } ##统计所有包含这个突变的sample(reads>=1)
	  
	  next if($#sample==0); ##过滤掉specific的site
	  
	  for $sample(@sample){
	  	$hash{$group{$sample}}++;
	  }
	  
	  
	  my @count=();my %temp=();
	  for $key(sort keys %hash){
	  	push(@count,$key);
	  }
	  @count=grep{++$temp{$_}<2} @count;
	  next if($#count==7);##过滤掉所有branch都有的site
	  
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
	  }
	  
	  next LABEL if($branch2>1);
	  next LABEL if($branch5>2);
	  next LABEL if($branch6>3);##tolerate 90% samples in every branch discard mutation(branch2:9 branch5:15 branch6:26)
	  
	  for $key(sort keys %hash){
	  	print INPUT "$key:$hash{$key}/$number{$key};";
	  }
	  print INPUT "\t$_\t@count\n";
	}
	close(DATA);
	close(INPUT);
}


@file=glob("D:/file/willow/newref/snp_newgroup/*.stats") or die "can't open glob";

open(INPUT,">D:/file/willow/newref/snp_newgroup/all.txt");
for $filename(@file){
	open(DATA,"<$filename");
	while(<DATA>){
		chomp;
		print INPUT "$_\n";
	}
	close(DATA);
}
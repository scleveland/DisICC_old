#!/usr/bin/perl
    
my @files = glob("*.fasta");

foreach my $file(@files)
{
   #print "FileName With Path[$file]\n";

	 $new_file = $file;
     $new_file =~ s/\..*//; #strip extension
     $new_file =~ s/.*\///; #strip path
     $new_file =~ s/.*\[//;
     $new_file =~ s/\].*//;
 
	 print "FileName Without Path $new_file\n";

	open (MYFILE, ">$new_file.fasta");

open(MYINPUTFILE, "<$file");
while(<MYINPUTFILE>)
 {
	# Good practice to store $_ value because
 	# subsequent operations may change it.
 	my($line) = $_;
 	if ($line=~ m/length/)
 	{
		my @fasta_line = split('length', $line);
		print MYFILE $fasta_line[0],"\n";
 	}
 	else
 	{
		print MYFILE $line;
 	}
 }

close(MYFILE);
}
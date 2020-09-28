#!/usr/bin/perl


open(OUT,">do_all_new.sh");
open(FH,"<get_mag_cubes_MUSE.csv");
while($line=<FH>) {
    chop($line);
    if ($line !~ "#") {
	my @data=split(",",$line);
	$name=$data[0];
	$red=$data[1];
	$check_file="/disk-d/manga/carlos/ccf_products/MUSE/ccf_dataproducts/".$name.".Ha_NII.ccf.vel.fits";
	if (!(-e $check_file)) {
	    $call="tsp ./ccf_IFS.py /mnt/synology/bacterio/disk-c/sanchez/MUSE/AMUSING/analysis/all_files/".$name."/GAS.".$name.".cube.fits.gz ".$red." /disk-d/manga/carlos/ccf_products/MUSE/ccf_dataproducts/".$name.".Ha_NII. 0";
	    print OUT "$call\n";
	} else {
	    print "$name already ccf analyzed\n";
	}
    }
}
close(FH);
close(OUT);

exit;

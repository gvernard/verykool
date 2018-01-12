<?php 

$out = $argv[1];


$fh = fopen($out . "cosmosis_output.txt","r");
$fh2 = fopen($out . "corner.txt","w");
while( $line = fgets($fh) ){
    if( substr($line,0,1) == '#' ){
	continue;
    } else {
	$single = preg_replace('!\s+!',' ',trim($line));
	$dum = explode(" ",$single);
	fprintf($fh2," 1 %12.5f",end($dum));
	for($i=0;$i<count($dum)-1;$i++){
	    fprintf($fh2," %12.5f",$dum[$i]);
	}
	fprintf($fh2,"\n");
    }
}
fclose($fh);
fclose($fh2);




?>

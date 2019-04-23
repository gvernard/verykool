<?php
// If 3 arguments are given (path,run,lmodel or step) it determines the 4th (step or lmodel)
// If 2 arguments are given (path,run), it determines the step and the lmodel
$path   = $argv[1];
$run    = $argv[2];


if( count($argv) == 3 ){

    // Determine lmodel 
    $lmodel = determineLmodel($path,$run);
    // Determine step
    $step = determineStep($path,$run,$lmodel);

} else if( count($argv) == 4 ){

    if( $argv[3] == "pert" or $argv[3] == "smooth" ){
	$lmodel = $argv[3];
	// Determine step
	$step = determineStep($path,$run,$lmodel);
    } else {
	$step = $argv[3];
	// Determine lmodel
	$lmodel = determineLmodel($path,$run);
    }

} else {
    
    $lmodel = "problem";
    $step   = "problem";

}

echo $lmodel . " " . $step;
#end of program


function determineLmodel($path,$run){
    if( file_exists($path . $run . "output/pert_initial_output.json") ){
	$lmodel = "pert";
    } else {
	$lmodel = "smooth";
    }
    return $lmodel;
}

function determineStep($path,$run,$lmodel){
    if( file_exists($path . $run . "output/".$lmodel."_model.fits") ){
	$step = "dum";
    } else {
	$dir = $path . $run . "output/";
	$scanned_dir = array_diff(scandir($dir),array('..','.'));
	$dum_arr = preg_grep("/".$lmodel."_model.fits/i",$scanned_dir);
	$dum_arr = array_values($dum_arr);
	
	$steps = array();
	foreach($dum_arr as $dum){
	    $tmp = explode('_',$dum);
	    if( is_numeric($tmp[0]) ){
		$steps[] = $tmp[0];
	    }
	}
	sort($steps);
	$step = end($steps);
    }
    return $step;
}

?>

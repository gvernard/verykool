<?php

$prog = $argv[1];
$path = $argv[2];
$run  = $argv[3];


// select suitable program output
$exec = "python";
switch($prog){
    case "plot_all.py":
	$out = "all.png";
	break;
    case "plot_all_pert.py":
	$out = "all.png";
	break;
    case "plot_corner.py":
	$out = "corner.pdf";
	break;
    case "plot_pow_spec.py":
	$out = "ps.png";
	break;
}


// scan run output directory to find the steps
$dir = $path . $run . "output/";
$scanned_dir = array_diff(scandir($dir),array('..','.'));
$dum_arr = preg_grep("/_model.fits/i",$scanned_dir);
$steps = array();
foreach($dum_arr as $dum){
    $tmp = explode('_',$dum);
    if( is_numeric($tmp[0]) ){
	$steps[] = $tmp[0];
    }
}
sort($steps);


// create output dir if it does not exist
$dir_out = "movie_" . $run;
if( file_exists($dir_out) ){
    exec('rm -r ' . $dir_out);
}
mkdir($dir_out);


// Loop over steps and produce outputs
echo "Looping over " . count($steps) . " steps using \"" . $prog . "\" :" . "\n";
foreach($steps as $step){
    echo $step . "\n";
    exec(implode(" ",array($exec,$prog,$path,$run,$step)));
    rename($out,$dir_out . "/" . str_pad($step,4,"0",STR_PAD_LEFT) . "_" . $out);
}

exec("ffmpeg -r 1 -i $dir_out/%4d_$out -r 10 $dir_out/video.mp4");
?>

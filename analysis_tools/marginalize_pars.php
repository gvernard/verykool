<?php 
# This is a dum program:
#   takes as input a file with columns: weights,like,par1,par2,...
#   and then computes the average and standard deviation for each parameter column, reading line by line
ini_set('memory_limit', '8192M'); 

$path = $argv[1];
$run  = $argv[2];
$sigma_interval = $argv[3];

if( count($argv) == 5 ){
    $lmodel = $argv[4];
    $step = "";
    $out_path = $path . $run . "output/" . $lmodel;
} else if( count($argv) == 6 ){
    $lmodel = $argv[4];
    $step = $argv[5];
    $out_path = $path . $run . "output/" . $step . "_" . $lmodel;
} else {
    echo "Either 4 or 5 command line arguments required: path, run, signam interval, lmodel, <step>\n";
    echo count($argv) . " provided, exiting!!!\n";
}




$i_prob   = 1;
$i_weight = -1;
$i_par    = 2;
if( $sigma_interval == 3 ){
    $Xsigma = 0.9973;
} else if( $sigma_interval == 2 ){
    $Xsigma = 0.9545;
} else {
    $Xsigma = 0.6827;
}


# Read parameter names and order
$par_names = array();
$fh   = fopen($path . $run . "output/" . $lmodel . "_postdist.paramnames","r");
while( $line = fgets($fh) ){
    $line = trim(preg_replace('!\s+!',' ',$line));
    $dum  = explode(' ',$line);
    $par_names[] = $dum[0];
    echo end($par_names)."\n";
}
fclose($fh);


# Read master file
$fh   = fopen($out_path . "_postdist.txt","r");
$line = fgets($fh);
$line = trim(preg_replace('!\s+!',' ',$line));
$dum  = explode(' ',$line);
$cols = count($dum);
fclose($fh);
printf("Number of columns: %d\n",$cols);

$lines = 0;
$fh = fopen($out_path . "_postdist.txt","r");
while( $line = fgets($fh) ){
    $lines++;
}
fclose($fh);
printf("Number of lines: %d\n",$lines);


# Weight probability accordingly
if( $i_weight < 0 ){

} else {
    echo " using weights\n";
    $weights = getColumn($out_path . "_postdist.txt",$i_weight,$lines);
    for($i=0;$i<count($weights);$i++){
	$prob[$i] *= $weights[$i];
    }
}





# Loop over parameters
$fout = fopen("sigma_intervals.dat","w");
fprintf($fout,"%20s %20s %20s %20s\n","name","mean","-" . $sigma_interval . " sigma","+" . $sigma_interval . " sigma");
for($i_par=0;$i_par<$cols-2;$i_par++){

    # Get probability. Need to do this every time because of the sorting
    $prob = getColumn($out_path . "_postdist.txt",$i_prob,$lines);


    # Sort probability simultaneously with parameter
    $x = getColumn($out_path . "_postdist.txt",$i_par+2,$lines);
    //for($i=0;$i<count($x);$i++){
    //    $x[$i] = log($x[$i]/25.0);
    //}
    array_multisort($x,$prob);
    

    # Reduce the arrays by adding the identical parameter values (arrays must be sorted first)
    $new_x = array();
    $new_p = array();
    reduceArrays($x,$prob,$new_x,$new_p);
    $N = count($new_x);
    

    # convert from probability density to probability
    $pbin = array_fill(0,$N,0);
    $pcum = array_fill(0,$N,0);
    for($i=1;$i<$N;$i++){
	$pbin[$i-1] = ($new_x[$i] + $new_x[$i-1])/2.0;
	$pcum[$i-1] = ($new_p[$i] + $new_p[$i-1])*($new_x[$i] - $new_x[$i-1])/2.0;
    }
    
    
    # convert to cumulative probability
    for($i=1;$i<$N;$i++){
	$pcum[$i] = $pcum[$i-1] + $pcum[$i];
    }
    

    # normalize
    for($i=0;$i<$N;$i++){
	$pcum[$i] /= $pcum[$N-1];
    }


    # Calculate mean, low and high interval
    $val  = 0.0;
    $low  = 0.0;
    $high = 0.0;
    
    $low_limit  = 0.5 - $Xsigma/2.0;
    $high_limit = 0.5 + $Xsigma/2.0;
    
    $val = calculateInterval($pbin,$pcum,0.5);
    if( $pcum[0] > $low_limit ){ // only an upper limit can be computed
	$high = calculateInterval($pbin,$pcum,$Xsigma);
    } else if( $pcum[$N-1] < $high_limit ) {
	$low = calculateInterval($pbin,$pcum,1.0 - $Xsigma);
    } else {
	$low  = calculateInterval($pbin,$pcum,$low_limit);
	$high = calculateInterval($pbin,$pcum,$high_limit);
    }


    fprintf($fout,"%20s %20.8f %20.8f %20.8f\n",$par_names[$i_par],$val,$low,$high);
}
fclose($fout);











function getColumn($fname,$col_index,$nlines){
    $vals = array_fill(0,$nlines,0);
    $fh = fopen($fname,"r");
    $i = 0;
    while( $line = fgets($fh) ){
	$line     = trim(preg_replace('!\s+!',' ',$line));
	$dum      = explode(' ',$line);
	$vals[$i] = floatval($dum[$col_index]);
	$i++;
    }
    fclose($fh);
    return $vals;
}

function reduceArrays($in_x,$in_p,&$out_x,&$out_p){
    $i = 0;
    $dum_x = array();
    $dum_p = array();
    while( $i<count($in_x) ){
	$dum_x[] = $in_x[$i];
	$val     = $in_p[$i];
	$j = $i+1;
	while( $j < count($in_x) and $in_x[$i] == $in_x[$j] ){
	    $val += $in_p[$j];
	    $j++;
	}
	$dum_p[] = $val;
	$i = $j;
    }
    $out_x = $dum_x;
    $out_p = $dum_p;
}

function calculateInterval($x,$y,$val){
    $x0;
    
    for($i=0;$i<count($y);$i++){
	if( $y[$i] >= $val ){
	    $x0 = $x[$i-1] + ($x[$i] - $x[$i-1])*($val - $y[$i-1])/($y[$i] - $y[$i-1]);
	    break;
	}
    }

    return $x0;
}




/*
$means = array_fill(0,$cols,0);
$prob = array();
$fh = fopen($fin,"r");
while( $line = fgets($fh) ){
    $line   = trim(preg_replace('!\s+!',' ',$line));
    $dum    = explode(' ',$line);
    $weight = floatval($dum[0]);
    $like   = floatval($dum[1]);
    $prob   = $like * $weight;
    $sum    += $prob;

    for($i=0;$i<$cols;$i++){
	$means[$i] += $prob * floatval($dum[$i+2]);
    }
}
fclose($fh);

for($i=0;$i<$cols;$i++){
    $means[$i] /= $sum;
}




$sdevs = array_fill(0,$cols,0);
$fh = fopen($fin,"r");
while( $line = fgets($fh) ){
    $line   = trim(preg_replace('!\s+!',' ',$line));
    $dum    = explode(' ',$line);
    $weight = floatval($dum[0]);
    $like   = floatval($dum[1]);
    $prob   = $like * $weight;

    for($i=0;$i<$cols;$i++){
	$sdevs[$i] += $prob * pow(floatval($dum[$i+2]) - $means[$i],2);
    }
}
fclose($fh);


printf("%12s %12s\n","means","sdevs");
for($i=0;$i<$cols;$i++){
    printf("%12.5e %12.5e\n",$means[$i],sqrt($sdevs[$i]/$sum));
}
*/

?>

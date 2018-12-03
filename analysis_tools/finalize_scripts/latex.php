<?php
$path = trim($argv[1]);
$run  = trim($argv[2]);

if( count($argv) == 4 ){
    $lmodel = $argv[3];
    $step = "";
    $out_path = $path . $run . "output/" . $lmodel;
} else if( count($argv) == 5 ){
    $lmodel = $argv[3];
    $step = $argv[4];
    $out_path = $path . $run . "output/" . $step . "_" . $lmodel;
} else {
    echo "Either 3 or 4 command line arguments required: path, run, lmodel, <step>\n";
    echo count($argv) . " provided, exiting!!!\n";
}



$dum = file_get_contents($path . $run .  'vkl_input.json');
$dum = preg_replace("#(/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/)|([\s\t]//.*)|(^//.*)#", '', $dum);
$inp = json_decode($dum,true);

$dum    = file_get_contents($out_path . '_output.json');
$out    = json_decode($dum,true);





$params = array();
foreach($inp['lenses'] as $key=>$lens){
    foreach($lens['nlpars'] as $param){
	if( $param['fix'] == 0 ){
	    //$param['map'] = $out['pars'][$param['nam']];
	    $params[$key . '_' . $param['nam']] = $param;
	}
    }
}
foreach($inp['reg']['nlpars'] as $param){
    if( $param['fix'] == 0 ){
	//$param['map'] = $out['pars'][$param['nam']];
	$params[$param['nam']] = $param;
    }
}
foreach($inp['physical']['nlpars'] as $param){
    if( $param['fix'] == 0 ){
	//$param['map'] = $out['pars'][$param['nam']];
	$params[$param['nam']] = $param;
    }
}


/*
$fh = fopen($path . 'output/vkl_stats.txt','r');
while( $dum = fscanf($fh,"%s%f%f%f%f\n") ){
    list($name,$mean,$sdev,$best,$map) = $dum;
    $params[$name]['mean'] = $mean;
    $params[$name]['sdev'] = $sdev;
    $params[$name]['best'] = $best;
    $params[$name]['map']  = $map;
}
fclose($fh);
*/




// Title
$str = '';
$str .= '\section*{\hfil Name: '.str_replace('_',' ',$run).' \hfil}';
file_put_contents('title.tex',$str);




// True, mean, sdev and MAP parameters
$str =  '';
/*
$str .= '\begin{table*}[!htb]';
$str .= '\caption{True vs Maximum A Posteriori (MAP) parameters.}';
$str .= '\begin{center}';
$str .= '\begin{tabular}{ r r r r r r r r }';
$str .= " " . implode(' & ',array('parameters','true','MAP','mean','sdev','min','max')) .  " \\\\ \n";
$str .= '\hline';

foreach($params as $param){
    $str .= " " . implode(' & ',array($param['nam'],$param['val'],round($param['map'],3),round($param['mean'],3),round($param['sdev'],3),$param['min'],$param['max'])) . " \\\\ \n";
}

$str .= '\hline';
$str .= '\end{tabular}';
$str .= '\end{center}';
$str .= '\label{tab:1}';
$str .= '\end{table*}';
*/
file_put_contents('table1.tex',$str);






// Table with the terms that contribute to the evidence calculation
$str =  '';
/*
$str .= '\begin{table*}[!htb]';
$str .= '\caption{Values of the various evidence terms at the MAP parameter values.}';
$str .= '\begin{center}';
$str .= '\begin{tabular}{ r l }';
$str .= ' term & value \\\\' . "\n";
$str .= '\hline';

$str .= ' $-E_D$                                                & ' . round(floatval(- $out["terms"]["ED"]),3) . ' \\\\ ' . "\n";
$str .= ' $-E_S$                                                & ' . round(floatval(- $out["terms"]["ES"]),3) . ' \\\\ ' . "\n";
$str .= ' $-\lambda E_S$                                        & ' . round(floatval($out["pars"]["lambda"])*floatval(- $out["terms"]["ES"]),3) . ' \\\\ ' . "\n";
$str .= ' $-E_D - \lambda E_S$                                  & ' . round(floatval(- $out["terms"]["ED"]) + floatval($out["pars"]["lambda"])*floatval(- $out["terms"]["ES"]),3) . ' \\\\ ' . "\n";
$str .= ' $-\frac{1}{2} \textrm{log} (\textrm{det} A )$         & ' . round(floatval($out["terms"]["logdetA"]),3) . ' \\\\ ' . "\n";
$str .= ' $ \frac{1}{2} \textrm{log} (\textrm{det} H^TH )$      & ' . round(floatval($out["terms"]["logdetHtH"]),3) . ' \\\\ ' . "\n";
$str .= ' $ \frac{1}{2} \textrm{log} (\textrm{det} C^{-1}_D )$  & ' . round(floatval($out["terms"]["logdetC"]),3) . ' \\\\ ' . "\n";
$str .= ' $ \frac{N_{source}}{2} \textrm{log} \lambda$          & ' . round(floatval($out["terms"]["logl"]),3) . ' \\\\ ' . "\n";
$str .= ' $-\frac{N_{data}}{2} \textrm{log} 2\pi$               & ' . round(floatval($out["terms"]["log2pi"]),3) . ' \\\\ ' . "\n";

$str .= '\hline';
$str .= '\end{tabular}';
$str .= '\end{center}';
$str .= '\label{tab:2}';
$str .= '\end{table*}';
*/
file_put_contents('table2.tex',$str);




// Table of other values of interest
if( file_exists($path . $run . "output/" . $lmodel . "_minimizer_output.json") ){
    $min = json_decode($path . $run . "output/" . $lmodel . "_minimizer_output.json",true);
    $iters = $min;
} else {
    $iters = array("n/a","n/a");
}


$str =  '';
$str .= '\begin{table*}[!htb]';
$str .= '\caption{Other values of interest.}';
$str .= '\begin{center}';
$str .= '\begin{tabular}{ r l l }';

$str .= '\hline';

$str .= ' $N_{source}$ & number of source pixels           & ' . $out["generic"]["Nsource"] . ' \\\\ ' . "\n";
$str .= ' $N_{data}$   & number of data pixels             & ' . $out["generic"]["Ndata"] . ' \\\\ ' . "\n";
if( strlen($inp['psfpath']) < 2 ){
    $str .= ' PSF & is PSF included?                          & No  \\\\' . "\n";
} else {
    $str .= ' PSF & is PSF included?                          & Yes \\\\' . "\n";
}
if( $inp['noise_flag'] == 'uniform' ){
    $cov = trim(file_get_contents($path . "data/noise.dat"));
    $str .= ' $C_D$        & covariance of the data            & uniform(' . $cov . ') \\\\ ' . "\n";
} else if( $inp['noise_flag'] == 'map' ){
    $str .= ' $C_D$        & covariance of the data            & map \\\\ ' . "\n";
}
$str .= ' $P_{data}$   & size of data pixels               & ' . sprintf("%10.5f",$out["generic"]["Psize"]) . ' \\\\ ' . "\n";
$str .= ' $N_{iter}$   & number of MultiNest iterations    & ' . $iters['total_samples'] . ' \\\\ ' . "\n";
$str .= ' $N_{repl}$   & number of live point replacements & ' . $iters['replacements'] . ' \\\\ ' . "\n";

$str .= '\hline';
$str .= '\end{tabular}';
$str .= '\end{center}';
$str .= '\label{tab:3}';
$str .= '\end{table*}';
file_put_contents('table3.tex',$str);





?>

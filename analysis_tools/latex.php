<?php
$path = trim($argv[1]);
$step = trim($argv[2]);

$dum   = explode('/',$path);
$case_name = $dum[count($dum)-2];

$dum    = file_get_contents($path . 'vkl_input.json');
$dum    = preg_replace("#(/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/)|([\s\t]//.*)|(^//.*)#", '', $dum);
$inp    = json_decode($dum,true);

$dum    = file_get_contents($path . 'output/'.$step.'_vkl_output.json');
$out    = json_decode($dum,true);





$params = array();
foreach($inp['lenses'][0]['nlpars'] as $param){
    if( $param['fix'] == 0 ){
	//$param['map'] = $out['pars'][$param['nam']];
	$params[$param['nam']] = $param;
    }
}
foreach($inp['reg']['nlpars'] as $param){
    if( $param['fix'] == 0 ){
	//$param['map'] = $out['pars'][$param['nam']];
	$params[$param['nam']] = $param;
    }
}
foreach($inp['physical'] as $param){
    if( $param['fix'] == 0 ){
	//$param['map'] = $out['pars'][$param['nam']];
	$params[$param['nam']] = $param;
    }
}


$fh = fopen($path . 'output/vkl_stats.txt','r');
while( $dum = fscanf($fh,"%s%f%f%f%f\n") ){
    list($name,$mean,$sdev,$best,$map) = $dum;
    $params[$name]['mean'] = $mean;
    $params[$name]['sdev'] = $sdev;
    $params[$name]['best'] = $best;
    $params[$name]['map']  = $map;
}
fclose($fh);





// Title
$str = '';
$str .= '\section*{\hfil Name: '.str_replace('_',' ',$case_name).' \hfil}';
file_put_contents('title.tex',$str);




// True, mean, sdev and MAP parameters
$str =  '';
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
file_put_contents('table1.tex',$str);






// Table with the terms that contribute to the evidence calculation
$str =  '';
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
file_put_contents('table2.tex',$str);




// Table of other values of interest
$fh    = fopen($path . "output/mn-resume.dat","r");
$line  = fgets($fh);
$line  = fgets($fh);
$line  = trim(preg_replace('!\s+!',' ',$line));
$iters = explode(' ',$line);
fclose($fh);



$str =  '';
$str .= '\begin{table*}[!htb]';
$str .= '\caption{Other values of interest.}';
$str .= '\begin{center}';
$str .= '\begin{tabular}{ r l l }';

$str .= '\hline';

$str .= ' $N_{source}$ & number of source pixels           & ' . $out["other"]["Nsource"] . ' \\\\ ' . "\n";
$str .= ' $N_{data}$   & number of data pixels             & ' . $out["other"]["Ndata"] . ' \\\\ ' . "\n";
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
$str .= ' $P_{data}$   & size of data pixels               & ' . $out["other"]["Psize"] . ' \\\\ ' . "\n";
$str .= ' $N_{iter}$   & number of MultiNest iterations    & ' . $iters[1] . ' \\\\ ' . "\n";
$str .= ' $N_{repl}$   & number of live point replacements & ' . $iters[0] . ' \\\\ ' . "\n";

$str .= '\hline';
$str .= '\end{tabular}';
$str .= '\end{center}';
$str .= '\label{tab:3}';
$str .= '\end{table*}';
file_put_contents('table3.tex',$str);





?>
